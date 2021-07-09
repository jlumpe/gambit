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

from gambit.io.seq import SequenceFile, find_kmers_in_files
from gambit.signatures.hdf5 import HDF5Signatures
from gambit.db.gambitdb import GAMBITDatabase
from gambit.db.models import ReferenceGenomeSet
from gambit.query import runquery
from gambit.cli import cli


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
		# Raises FileNotFoundError if does not exist
		files[k] = v.resolve()

	return files


@pytest.fixture(scope='module')
def signatures(testdb_files):
	"""K-mer signatures for test genomes."""
	return HDF5Signatures.open(testdb_files['signatures'])


@pytest.fixture(scope='module')
def testdb(testdb_session, signatures):
	"""Full GAMBITDatabase object for test db."""

	session = testdb_session()
	gset = session.query(ReferenceGenomeSet).one()

	return GAMBITDatabase(gset, signatures)


@pytest.fixture(scope='module')
def query_data(testdb_files):
	"""Query files and their expected taxa."""
	table_path = testdb_files['queries']
	seqs_dir = table_path.parent

	files = []
	expected_taxa = []

	with open(table_path, newline='') as f:
		for row in DictReader(f):
			file = SequenceFile(
				path=seqs_dir / (row['name'] + '.fa'),
				format='fasta',
			)
			files.append(file)
			expected_taxa.append(row['expected_taxon'])

	return files, expected_taxa


def test_query_python(testdb, query_data):
	"""Run a full query using the Python API."""

	query_files, expected_taxa = query_data
	query_sigs = find_kmers_in_files(testdb.kmerspec, query_files)

	results = runquery(testdb, query_sigs, query_files)

	assert len(results.items) == len(query_files)

	for item, file, expected_key in zip(results.items, query_files, expected_taxa):
		expected_taxon = testdb.genomeset.taxa.filter_by(key=expected_key).one() if expected_key else None

		assert item.input.file == file
		assert item.success

		if expected_taxon is None:
			assert item.predicted_taxon is None
			assert item.report_taxon is None
			assert item.primary_match is None
		else:
			assert item.predicted_taxon == expected_taxon
			assert item.report_taxon == expected_taxon
			assert item.primary_match is not None
			assert item.primary_match.matching_taxon == expected_taxon

			# In this database, closest match should be primary match
			assert item.closest_match == item.primary_match

		assert not item.warnings
		assert item.error is None


@pytest.mark.parametrize('out_fmt', ['json'])
def test_query_cli(testdb_files, query_data, out_fmt, tmp_path):
	"""Run a full query using the command line interface."""
	results_file = tmp_path / 'results.json'
	query_files, expected_taxa = query_data

	args = [
		f'--db={testdb_files["root"]}',
		'query',
		f'--output={results_file}',
		f'--outfmt={out_fmt}',
		*(str(f.path) for f in query_files),
	]

	cli.main(args, standalone_mode=False)

	if out_fmt == 'json':
		_check_results_json(results_file, query_files, expected_taxa)
	else:
		assert False

def _check_results_json(results_file, query_files, expected_taxa):
	with results_file.open() as f:
		results = json.load(f)

	items = results['items']
	assert isinstance(items, list)
	assert len(items) == len(query_files)

	for item, file, expected in zip(items, query_files, expected_taxa):
		assert item['input']['label'] == file.path.name
		assert item['success'] is True

		if expected == '':
			assert item['predicted_taxon'] is None
			assert item['report_taxon'] is None
			assert item['primary_match'] is None
		else:
			predicted = item['predicted_taxon']
			assert predicted is not None
			assert predicted['key'] == expected
			assert item['report_taxon'] == predicted
			assert item['primary_match'] is not None
			assert item['primary_match']['matching_taxon']['key'] == expected

			# In this database, closest match should be primary match
			assert item['closest_match'] == item['primary_match']

		assert item['warnings'] == []
		assert item['error'] is None
