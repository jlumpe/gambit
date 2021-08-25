"""
Test the 'gambit query' CLI command using the testdb_210818 database.
"""

import json
from csv import DictReader
from io import StringIO
from pathlib import Path

import pytest
import numpy as np

from gambit.cli import cli
from gambit.util.misc import zip_strict
from gambit.io.export.json import JSONResultsExporter
from gambit.io.export.csv import CSVResultsExporter
from gambit import __version__ as GAMBIT_VERSION


@pytest.mark.parametrize('out_fmt', ['csv', 'json'])
def test_query_cmd(testdb_files, testdb, testdb_queries, testdb_results, out_fmt, tmp_path):
	"""Run a full query using the command line interface."""
	results_file = tmp_path / ('results.' + out_fmt)
	query_files = [query['file'] for query in testdb_queries]
	params = testdb_results.params

	args = [
		f'--db={testdb_files["root"]}',
		'query',
		f'--output={results_file}',
		f'--outfmt={out_fmt}',
		'--strict' if params.classify_strict else '--no-strict',
		*(str(f.path) for f in query_files),
	]

	cli.main(args, standalone_mode=False)

	# Detailed checks of output format are already present in tests for exporter classes, just check
	# that the exported data matches an export of the reference results
	if out_fmt == 'json':
		_check_results_json(results_file, query_files, testdb_results)
	elif out_fmt == 'csv':
		_check_results_csv(results_file, query_files, testdb_results)
	else:
		assert False


def _check_results_json(results_file, query_files, ref_results):
	with results_file.open() as f:
		data = json.load(f)

	# Equivalent data for reference results
	exporter = JSONResultsExporter()
	buf = StringIO()
	exporter.export(buf, ref_results)
	buf.seek(0)
	ref_data = json.load(buf)

	assert data['gambit_version'] == GAMBIT_VERSION
	assert len(data['items']) == len(query_files)

	for key in ['genomeset', 'signaturesmeta', 'extra']:
		assert data[key] == ref_data[key]

	for item, ref_item, query_file in zip_strict(data['items'], ref_data['items'], query_files):
		assert Path(item['query']['path']).name == query_file.path.name
		assert item['query']['format'] == query_file.format

		assert item['predicted_taxon'] == ref_item['predicted_taxon']
		assert item['closest_genome'] == ref_item['closest_genome']
		assert np.isclose(item['closest_genome_distance'], ref_item['closest_genome_distance'])


def _check_results_csv(results_file, query_files, ref_results):
	with results_file.open() as f:
		rows = list(DictReader(f))

	exporter = CSVResultsExporter()
	buf = StringIO()
	exporter.export(buf, ref_results)
	buf.seek(0)
	ref_rows = list(DictReader(buf))

	assert len(rows) == len(ref_rows)

	cmp_cols = [
		'predicted.name',
		'predicted.rank',
		'predicted.ncbi_id',
		'predicted.threshold',
		'closest.description',
	]

	for row, ref_row, file in zip_strict(rows, ref_rows, query_files):
		assert row['query.path'] == str(file.path)
		assert np.isclose(float(row['closest.distance']), float(ref_row['closest.distance']))

		for key in cmp_cols:
			assert row[key] == ref_row[key]
