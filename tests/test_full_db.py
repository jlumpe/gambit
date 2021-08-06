"""
Run a full set of queries using the testdb_210126 database.

Tests in this file will only be run when the --gambit-test-full-db option is passed to the pycharm
command.

Database files and query sequences are located in tests/data/testdb_210126, but only the genome
database file is included in version control. Other files need to be obtained separately. See the
Readme.md file in that directory for more information.
"""

import json
from csv import DictReader

import pytest

from gambit.io.seq import SequenceFile
from gambit.signatures.hdf5 import HDF5Signatures
from gambit.db import GAMBITDatabase, ReferenceGenomeSet
from gambit.query import QueryParams, query_parse
from gambit.cli import cli
from gambit.util.misc import zip_strict
from gambit import __version__ as GAMBIT_VERSION


@pytest.fixture(autouse=True, scope='module')
def testdb_files(request, testdb_dir):
	"""Paths to testdb_210126 files.

	Skips all dependent tests if the --gambit-test-full-db command line option is not passed.

	Checks that the directory and required files/subdirectories exist and fails tests immediately
	if they do not.
	"""
	if not request.config.getoption('gambit_test_full_db'):
		pytest.skip('--gambit-test-full-db option not given')

	files = dict(
		root=testdb_dir,
		db=testdb_dir / 'testdb_210126-genomes.db',
		signatures=testdb_dir / 'testdb_210126-signatures.h5',
		queries=testdb_dir / 'query-seqs/queries.csv',
	)

	for k, v in files.items():
		assert v.exists(), f'Required testdb file not found: {v}'

	return files


@pytest.fixture(scope='module')
def signatures(testdb_files):
	"""K-mer signatures for test genomes."""
	return HDF5Signatures.open(testdb_files['signatures'])


@pytest.fixture(scope='module')
def testdb(testdb_session, signatures):
	"""Full GAMBITDatabase object for test db."""

	with testdb_session() as session:
		gset = session.query(ReferenceGenomeSet).one()
		yield GAMBITDatabase(gset, signatures)


@pytest.fixture(scope='module')
def query_data(testdb, testdb_files):
	"""Query files and their expected taxa."""
	table_path = testdb_files['queries']
	seqs_dir = table_path.parent

	files = []
	expected_taxa = []

	with open(table_path, newline='') as f:
		for row in DictReader(f):
			files.append(SequenceFile(
				path=seqs_dir / (row['name'] + '.fa'),
				format='fasta',
			))

			if row['expected_taxon']:
				expected = testdb.genomeset.taxa.filter_by(key=row['expected_taxon']).one()
			else:
				expected = None
			expected_taxa.append(expected)

	return files, expected_taxa


@pytest.mark.parametrize('classify_strict', [False, True])
def test_query_python(testdb, query_data, classify_strict):
	"""Run a full query using the Python API."""

	query_files, expected_taxa = query_data
	params = QueryParams(classify_strict=classify_strict)

	results = query_parse(testdb, query_files, params)

	assert results.params == params
	assert results.genomeset == testdb.genomeset
	assert results.signaturesmeta == testdb.signatures.meta
	assert results.gambit_version == GAMBIT_VERSION

	for item, file, expected_taxon in zip_strict(results.items, query_files, expected_taxa):
		clsresult = item.classifier_result

		assert item.input.file == file
		assert clsresult.success

		if expected_taxon is None:
			assert clsresult.predicted_taxon is None
			assert clsresult.primary_match is None
			assert item.report_taxon is None
		else:
			assert clsresult.predicted_taxon == expected_taxon
			assert item.report_taxon == expected_taxon
			assert clsresult.primary_match is not None
			assert clsresult.primary_match.matched_taxon == expected_taxon

			# In this database, closest match should be primary match
			assert clsresult.closest_match == clsresult.primary_match

		assert not clsresult.warnings
		assert clsresult.error is None


@pytest.mark.parametrize('out_fmt', ['csv', 'json'])
@pytest.mark.parametrize('classify_strict', [False, True])
def test_query_cli(testdb_files, testdb, query_data, out_fmt, classify_strict, tmp_path):
	"""Run a full query using the command line interface."""
	results_file = tmp_path / 'results.json'
	query_files, expected_taxa = query_data

	args = [
		f'--db={testdb_files["root"]}',
		'query',
		f'--output={results_file}',
		f'--outfmt={out_fmt}',
		'--strict' if classify_strict else '--no-strict',
		*(str(f.path) for f in query_files),
	]

	cli.main(args, standalone_mode=False)

	# Detailed checks of output format are already present in tests for exporter classes, just need
	# to check that the results themselves seem correct
	if out_fmt == 'json':
		_check_results_json(results_file, testdb, query_files, expected_taxa)
	elif out_fmt == 'csv':
		_check_results_csv(results_file, query_files, expected_taxa)
	else:
		assert False

def _check_results_json(results_file, testdb, query_files, expected_taxa):
	with results_file.open() as f:
		results = json.load(f)

	assert results['genomeset']['key'] == testdb.genomeset.key
	assert results['signaturesmeta']['id'] == testdb.signatures.meta.id
	assert results['gambit_version'] == GAMBIT_VERSION

	items = results['items']
	assert len(items) == len(query_files)

	for item, file, expected in zip_strict(items, query_files, expected_taxa):
		assert item['query']['path'] == str(file.path)

		if expected is None:
			assert item['predicted_taxon'] is None
		else:
			assert item['predicted_taxon']['key'] == expected.key

def _check_results_csv(results_file, query_files, expected_taxa):
	with results_file.open() as f:
		rows = list(DictReader(f))

	assert len(rows) == len(query_files)

	for row, file, expected in zip_strict(rows, query_files, expected_taxa):
		assert row['query.path'] == str(file.path)

		if expected is None:
			assert row['predicted.name'] == ''
		else:
			assert row['predicted.name'] == expected.name
