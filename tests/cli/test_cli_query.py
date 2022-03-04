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
	(Currently not overridden anywhere).

	Based on this example:
	https://docs.pytest.org/en/6.2.x/fixture.html#using-markers-to-pass-data-to-fixtures

	Note than with slice notation, `[:None]` is the same as `[:]`.
	"""
	marker = request.node.get_closest_marker("testdb_nqueries")
	return None if marker is None else marker.args[0]

@pytest.fixture()
def query_files(testdb_query_files, nqueries):
	"""Paths to query files."""
	return testdb_query_files[:nqueries]

@pytest.fixture(name='make_args')
def make_args_factory(testdb_files):

	def make_args(query_files=None, sigfile=None, output=None, outfmt=None, strict=False):
		"""Make command line arguments for query file."""

		args = [f'--db={testdb_files["root"]}', 'query']
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

	return make_args


def check_results(results_file, out_fmt, ref_results):
	"""Check results output matches reference QueryResults object."""
	if out_fmt == 'json':
		with open(results_file) as fh:
			check_json_results(fh, ref_results, strict=False)
	elif out_fmt == 'csv':
		with open(results_file) as fh:
			check_csv_results(fh, ref_results, strict=False)
	else:
		raise ValueError(f'Invalid out_fmt {out_fmt!r}')


@pytest.mark.parametrize(
	['out_fmt', 'testdb_results', 'testdb_queries_gzipped'],
	[
		('json', 'non_strict', False),
		('csv',  'non_strict', False),
		('json', 'strict',     False),
		('csv',  'strict',     False),
		('json', 'non_strict', True),
	],
	indirect=['testdb_results', 'testdb_queries_gzipped'],
)
def test_full_query(make_args, testdb_results, nqueries, query_files, out_fmt, tmp_path):
	"""Run a full query using the command line interface."""

	# Modify results object to match file names and possibly reduced # of queries
	ref_results = copy(testdb_results)
	ref_results.items = ref_results.items[:nqueries]
	for item, file in zip_strict(ref_results.items, query_files):
		item.input = QueryInput.convert(file)

	results_file = tmp_path / ('results.' + out_fmt)

	args = make_args(
		query_files=query_files,
		output=results_file,
		outfmt=out_fmt,
		strict=ref_results.params.classify_strict,
	)

	result = invoke_cli(args)
	assert result.exit_code == 0

	check_results(results_file, out_fmt, ref_results)


# Not really necessary to check all combinations of parameters.
@pytest.mark.parametrize('out_fmt', ['json'])
@pytest.mark.parametrize('testdb_results', ['non_strict'], indirect=True)
def test_sigfile(make_args, testdb_files, testdb_query_signatures, testdb_results, out_fmt, tmp_path):
	"""Test using signature file instead of parsing genome files."""

	results_file = tmp_path / ('results.' + out_fmt)

	# Modify results object inputs to match signature file
	ref_results = copy(testdb_results)
	for item, id_ in zip_strict(ref_results.items, testdb_query_signatures.ids):
		item.input = QueryInput(id_)

	args = make_args(
		sigfile=testdb_files['query_signatures'],
		output=results_file,
		outfmt=out_fmt,
		strict=ref_results.params.classify_strict,
	)

	result = invoke_cli(args)
	assert result.exit_code == 0

	check_results(results_file, out_fmt, ref_results)
