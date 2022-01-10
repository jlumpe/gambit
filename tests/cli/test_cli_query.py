"""
Test the 'gambit query' CLI command using the testdb_210818 database.
"""

from copy import copy

import pytest

from gambit.cli.test import invoke_cli
from gambit.results.test import check_json_results, check_csv_results
from gambit.query import QueryInput
from gambit.util.misc import zip_strict


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


def check_results(results_file, out_fmt,ref_results):
	"""Check results output matches reference QueryResults object."""
	if out_fmt == 'json':
		with open(results_file) as fh:
			check_json_results(fh, ref_results, strict=False)
	elif out_fmt == 'csv':
		with open(results_file) as fh:
			check_csv_results(fh, ref_results, strict=False)
	else:
		raise ValueError(f'Invalid out_fmt {out_fmt!r}')


@pytest.mark.parametrize('out_fmt', ['csv', 'json'])
def test_full_query(testdb_files, results, query_files, out_fmt, tmp_path):
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

	check_results(results_file, out_fmt, results)


@pytest.mark.parametrize('out_fmt', ['csv', 'json'])
def test_sigfile(testdb_files, testdb_query_signatures, results, out_fmt, tmp_path):
	"""Test using signature file instead of parsing genome files."""

	results_file = tmp_path / ('results.' + out_fmt)

	# Modify results object inputs to match signature file
	for item, id_ in zip_strict(results.items, testdb_query_signatures.ids):
		item.input = QueryInput(id_)

	args = make_args(
		sigfile=testdb_files['query_signatures'],
		db=testdb_files['root'],
		output=results_file,
		outfmt=out_fmt,
		strict=results.params.classify_strict,
	)

	result = invoke_cli(args)
	assert result.exit_code == 0

	check_results(results_file, out_fmt, results)
