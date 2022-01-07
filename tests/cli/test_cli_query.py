"""
Test the 'gambit query' CLI command using the testdb_210818 database.
"""

import json
from csv import DictReader
from copy import copy

import pytest
import numpy as np

from gambit.cli.test import invoke_cli
from gambit.io.results.base import export_to_buffer
from gambit.io.results.json import JSONResultsExporter
from gambit.io.results.csv import CSVResultsExporter


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


def make_args(query_files=None, sigfile=None, db=None, output=None, outfmt=None, strict=False):
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
	if sigfile is not None:
		args.append(f'--sigfile={sigfile}')

	if query_files is not None:
		args.extend([str(f.path) for f in query_files])

	return args


def check_results(results_file,
                  query_files,
                  out_fmt,
                  ref_results,
                  check_input_path=False,
                  input_labels=None,
                  ):
	"""Check results output matches reference QueryResults object.

	Detailed checks of output format are already present in tests for exporter classes, just check
	that the exported data matches an export of the reference results.
	"""
	if out_fmt == 'json':
		_check_results_json(results_file, query_files, ref_results, check_input_path=check_input_path, input_labels=input_labels)
	elif out_fmt == 'csv':
		_check_results_csv(results_file, query_files, ref_results, check_input_path=check_input_path, input_labels=input_labels)
	else:
		raise ValueError(f'Invalid out_fmt {out_fmt!r}')

def _check_results_json(results_file,
                        query_files,
                        ref_results,
                        check_input_path=False,
                        input_labels=None,
                        ):
	with results_file.open() as f:
		data = json.load(f)

	# Equivalent data for reference results
	buf = export_to_buffer(ref_results, JSONResultsExporter())
	ref_data = json.load(buf)

	assert len(data['items']) == len(query_files)

	for key in ['genomeset', 'signaturesmeta', 'extra']:
		assert data[key] == ref_data[key]

	# for item, ref_item, query_file in zip_strict(data['items'], ref_data['items'], query_files):
	for i, item in enumerate(data['items']):
		if check_input_path:
			assert item['query']['path'] == str(query_files[i].path)
		if input_labels is not None:
			assert item['query']['name'] == input_labels[i]
		assert item['query']['format'] == query_files[i].format

		ref_item = ref_data['items'][i]
		assert item['predicted_taxon'] == ref_item['predicted_taxon']
		assert item['closest_genome'] == ref_item['closest_genome']
		assert np.isclose(item['closest_genome_distance'], ref_item['closest_genome_distance'])

def _check_results_csv(results_file,
                       query_files,
                       ref_results,
                       check_input_path=False,
                       input_labels=None,
                       ):
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

	# for row, ref_row, file in zip_strict(rows, ref_rows, query_files):
	for i, row in enumerate(rows):
		if check_input_path:
			assert row['query.path'] == str(query_files[i].path)
		if input_labels is not None:
			assert row['query.name'] == input_labels[i]

		assert np.isclose(float(row['closest.distance']), float(ref_rows[i]['closest.distance']))

		for key in cmp_cols:
			assert row[key] == ref_rows[i][key]


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
		query_files=query_files,
		db=testdb_files['root'],
		output=results_file,
		outfmt=out_fmt,
		strict=results.params.classify_strict,
	)

	result = invoke_cli(args)
	assert result.exit_code == 0

	check_results(results_file, query_files, out_fmt, results)


def test_sigfile(testdb_files,
                 query_files,
                 testdb_query_signatures,
                 results,
                 tmp_path,
                 ):
	"""Test using signature file instead of parsing genome files."""

	out_fmt = 'csv'
	results_file = tmp_path / ('results.' + out_fmt)

	args = make_args(
		sigfile=testdb_files['query_signatures'],
		db=testdb_files['root'],
		output=results_file,
		outfmt=out_fmt,
		strict=results.params.classify_strict,
	)

	result = invoke_cli(args)
	assert result.exit_code == 0

	check_results(
		results_file,
		query_files,
		out_fmt,
		results,
		check_input_path=False,
		input_labels=testdb_query_signatures.ids,
	)
