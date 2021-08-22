"""
Test the 'gambit query' CLI command using the testdb_210818 database.
"""

import json
from csv import DictReader

import pytest

from gambit.cli import cli
from gambit.util.misc import zip_strict
from gambit import __version__ as GAMBIT_VERSION


@pytest.mark.parametrize('out_fmt', ['csv', 'json'])
@pytest.mark.parametrize('classify_strict', [False, True])
def test_query_cmd(testdb_files, testdb, testdb_queries, out_fmt, classify_strict, tmp_path):
	"""Run a full query using the command line interface."""
	if not classify_strict:
		pytest.skip()  # TODO

	results_file = tmp_path / 'results.json'
	query_files = [query['file'] for query in testdb_queries]

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
		_check_results_json(results_file, testdb, testdb_queries)
	elif out_fmt == 'csv':
		_check_results_csv(results_file, testdb_queries)
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
