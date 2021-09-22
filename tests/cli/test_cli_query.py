"""
Test the 'gambit query' CLI command using the testdb_210818 database.
"""

import json
from csv import DictReader
from io import StringIO
from pathlib import Path
from copy import copy

import pytest
import numpy as np
from click.testing import CliRunner

from gambit.cli import cli
from gambit.util.misc import zip_strict
from gambit.io.export.json import JSONResultsExporter
from gambit.io.export.csv import CSVResultsExporter


@pytest.fixture()
def nqueries(request):
	"""Number of testdb query files to use, None means use all of them.

	Value is derived from argument to the "testdb_nqueries" marker, if any. This can be set on
	specific test functions to improve speed by only using a subset of the query files.

	Based on this example:
	https://docs.pytest.org/en/6.2.x/fixture.html#using-markers-to-pass-data-to-fixtures
	"""
	marker = request.node.get_closest_marker("testdb_nqueries")
	return None if marker is None else marker.args[0]

@pytest.fixture()
def query_files(testdb_queries, nqueries):
	"""Paths to query files."""
	files = [query['file'] for query in testdb_queries]
	return files if nqueries is None else files[:nqueries]

@pytest.fixture()
def results(testdb_results, nqueries):
	"""Results object to compare output to."""
	if nqueries is None:
		return testdb_results

	results = copy(testdb_results)
	results.items = results.items[:nqueries]
	return results


def make_args(query_files, db=None, output=None, outfmt=None, strict=False):
	"""Make command line arguments for query file."""
	args = []

	if db is not None:
		args.append(f'--db={db}')

	args.append('query')
	args.append('--strict' if strict else '--no-strict')

	if output is not None:
		args.append(f'--output={output}')
	if outfmt is not None:
		args.append(f'--outfmt={outfmt}')

	args.extend([str(f.path) for f in query_files])

	return args


def export_to_buffer(results, exporter):
	buf = StringIO()
	exporter.export(buf, results)
	buf.seek(0)
	return buf


def check_results(results_file, query_files, out_fmt, ref_results):
	"""Check results output matches reference QueryResults object.

	Detailed checks of output format are already present in tests for exporter classes, just check
	that the exported data matches an export of the reference results.
	"""
	if out_fmt == 'json':
		_check_results_json(results_file, query_files, ref_results)
	elif out_fmt == 'csv':
		_check_results_csv(results_file, query_files, ref_results)
	else:
		raise ValueError(f'Invalid out_fmt {out_fmt!r}')

def _check_results_json(results_file, query_files, ref_results):
	with results_file.open() as f:
		data = json.load(f)

	# Equivalent data for reference results
	buf = export_to_buffer(ref_results, JSONResultsExporter())
	ref_data = json.load(buf)

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

	buf = export_to_buffer(ref_results, CSVResultsExporter())
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


@pytest.mark.parametrize('out_fmt', ['csv', 'json'])
def test_full_query(testdb_files,
                    query_files,
                    results,
                    out_fmt,
                    tmp_path,
                    ):
	"""Run a full query using the command line interface."""

	results_file = tmp_path / ('results.' + out_fmt)

	args = make_args(
		query_files,
		db=testdb_files["root"],
		output=results_file,
		outfmt=out_fmt,
		strict=results.params.classify_strict,
	)

	runner = CliRunner()
	result = runner.invoke(cli, args)
	assert result.exit_code == 0

	check_results(results_file, query_files, out_fmt, results)


@pytest.mark.testdb_nqueries(10)
def test_db_from_env(testdb_files,
                     query_files,
                     results,
                     tmp_path,
                     ):
	"""Test setting the database directory using an environment variable instead of an argument."""

	out_fmt = 'csv'
	results_file = tmp_path / ('results.' + out_fmt)

	args = make_args(
		query_files,
		output=results_file,
		outfmt=out_fmt,
		strict=results.params.classify_strict,
	)

	env = dict(
		GAMBIT_DB_PATH=str(testdb_files["root"]),
	)

	runner = CliRunner()
	result = runner.invoke(cli, args, env=env)
	assert result.exit_code == 0

	check_results(results_file, query_files, out_fmt, results)
