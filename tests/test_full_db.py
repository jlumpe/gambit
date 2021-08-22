"""
Run a full set of queries using the testdb_210818 database.

Database files and query sequences are located in tests/data/testdb_210818.
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


@pytest.fixture(scope='module')
def signatures(testdb_files):
	"""K-mer signatures for test genomes."""
	return HDF5Signatures.open(testdb_files['ref_signatures'])


@pytest.fixture(scope='module')
def testdb(testdb_session, signatures):
	"""Full GAMBITDatabase object for test db."""

	with testdb_session() as session:
		gset = session.query(ReferenceGenomeSet).one()
		yield GAMBITDatabase(gset, signatures)


@pytest.fixture(scope='module')
def queries(testdb, testdb_files):
	"""Query files and their expected results."""
	genomes_dir = testdb_files['query_genomes']

	with open(testdb_files['queries_table'], newline='') as f:
		rows = list(DictReader(f))

	for row in rows:
		row['warnings'] = row['warnings'].lower() == 'true'
		row['file'] = SequenceFile(
			path=genomes_dir / (row['name'] + '.fasta'),
			format='fasta',
		)

	return rows


@pytest.mark.parametrize('classify_strict', [False, True])
def test_query_python(testdb, queries, classify_strict):
	"""Run a full query using the Python API."""

	query_files = [item['file'] for item in queries]
	params = QueryParams(classify_strict=classify_strict)

	results = query_parse(testdb, query_files, params)

	assert results.params == params
	assert results.genomeset == testdb.genomeset
	assert results.signaturesmeta == testdb.signatures.meta
	assert results.gambit_version == GAMBIT_VERSION

	for query, item in zip_strict(queries, results.items):
		clsresult = item.classifier_result
		predicted = clsresult.predicted_taxon

		assert item.input.file == query['file']
		assert clsresult.success
		assert clsresult.error is None

		if classify_strict:
			if query['predicted']:
				assert predicted is not None
				assert predicted.name == query['predicted']
				assert clsresult.primary_match is not None
				assert clsresult.primary_match.genome.description == query['primary']
				assert item.report_taxon is (predicted if predicted.report else predicted.parent)

			else:
				assert predicted is None
				assert clsresult.primary_match is None
				assert item.report_taxon is None

			assert clsresult.closest_match.genome.description == query['closest']
			assert bool(clsresult.warnings) == query['warnings']

		else:
			if query['predicted']:
				assert clsresult.primary_match == clsresult.closest_match
				assert predicted is clsresult.primary_match.matched_taxon
				assert item.report_taxon is (predicted if predicted.report else predicted.parent)

			else:
				assert predicted is None
				assert clsresult.primary_match is None
				assert item.report_taxon is None

			assert clsresult.closest_match.genome.description == query['closest']
			assert not clsresult.warnings


@pytest.mark.parametrize('out_fmt', ['csv', 'json'])
@pytest.mark.parametrize('classify_strict', [False, True])
def test_query_cli(testdb_files, testdb, queries, out_fmt, classify_strict, tmp_path):
	"""Run a full query using the command line interface."""
	if not classify_strict:
		pytest.skip()  # TODO

	results_file = tmp_path / 'results.json'
	query_files = [query['file'] for query in queries]

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
		_check_results_json(results_file, testdb, queries)
	elif out_fmt == 'csv':
		_check_results_csv(results_file, queries)
	else:
		assert False

def _check_results_json(results_file, testdb, queries):
	with results_file.open() as f:
		results = json.load(f)

	assert results['genomeset']['key'] == testdb.genomeset.key
	assert results['signaturesmeta']['id'] == testdb.signatures.meta.id
	assert results['gambit_version'] == GAMBIT_VERSION

	items = results['items']
	assert len(items) == len(queries)

	for item, query in zip_strict(items, queries):
		assert item['query']['path'] == str(query['file'].path)

		if query['predicted']:
			assert item['predicted_taxon']['name'] == query['predicted']
		else:
			assert item['predicted_taxon'] is None

def _check_results_csv(results_file, queries):
	with results_file.open() as f:
		rows = list(DictReader(f))

	assert len(rows) == len(queries)

	for row, query in zip_strict(rows, queries):
		assert row['query.path'] == str(query['file'].path)

		if query['predicted']:
			assert row['predicted.name'] == query['predicted']
		else:
			assert row['predicted.name'] == ''
