"""
Test the 'gambit query' CLI command using the testdb_210818 database.
"""

from copy import copy

import pytest

from gambit.cli.test import invoke_cli
from gambit.results.test import check_json_results, check_csv_results
from gambit.query import QueryInput
from gambit.util.misc import zip_strict
from gambit.util.io import write_lines


@pytest.fixture(params=[None])
def nqueries(request):
	"""Number of testdb query files to use, None means use all of them.

	Can be changed via indirect parameterization in specific tests.
	Note than with slice notation, `[:None]` is the same as `[:]`.
	"""
	return request.param

@pytest.fixture()
def query_files(testdb_query_files, nqueries):
	"""Paths to query files."""
	return testdb_query_files[:nqueries]

@pytest.fixture(name='make_args')
def make_args_factory(testdb_files, query_files, tmp_path):

	def make_args(positional=False, list_file=False, sig_file=False, output=None, outfmt=None, strict=False):
		"""Make command line arguments for query file."""

		args = [f'--db={testdb_files["root"]}', 'query']
		args.append('--strict' if strict else '--no-strict')

		if output is not None:
			args.append(f'--output={output}')
		if outfmt is not None:
			args.append(f'--outfmt={outfmt}')

		if positional:
			args.extend([str(f.path) for f in query_files])
		if list_file:
			list_file = tmp_path / 'genomes.txt'
			write_lines([f.path.name for f in query_files], list_file)
			args += ['-l', str(list_file), f'--ldir={query_files[0].path.parent}']
		if sig_file:
			args.append(f'--sigfile={testdb_files["query_signatures"]}')

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
	['nqueries', 'use_list_file', 'out_fmt', 'testdb_results', 'testdb_queries_gzipped'],
	[
		(None, False, 'json', 'non_strict', False),
		(20,   False, 'csv',  'non_strict', False),
		(None, False, 'json', 'strict',     False),
		(20,   False, 'csv',  'strict',     False),
		(None, False, 'json', 'non_strict', True),
		(20,   True,  'json', 'non_strict', False),
	],
	indirect=['nqueries', 'testdb_results', 'testdb_queries_gzipped'],
)
def test_full_query(make_args, testdb_results, nqueries, use_list_file, query_files, out_fmt, tmp_path):
	"""Run a full query using the command line interface."""

	# Modify results object to match file names and possibly reduced # of queries
	ref_results = copy(testdb_results)
	ref_results.items = ref_results.items[:nqueries]
	for item, file in zip_strict(ref_results.items, query_files):
		item.input = QueryInput.convert(file)

	results_file = tmp_path / ('results.' + out_fmt)

	args = make_args(
		positional=not use_list_file,
		list_file=use_list_file,
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
def test_sigfile(make_args, testdb_query_signatures, testdb_results, out_fmt, tmp_path):
	"""Test using signature file instead of parsing genome files."""

	results_file = tmp_path / ('results.' + out_fmt)

	# Modify results object inputs to match signature file
	ref_results = copy(testdb_results)
	for item, id_ in zip_strict(ref_results.items, testdb_query_signatures.ids):
		item.input = QueryInput(id_)

	args = make_args(
		sig_file=True,
		output=results_file,
		outfmt=out_fmt,
		strict=ref_results.params.classify_strict,
	)

	result = invoke_cli(args)
	assert result.exit_code == 0

	check_results(results_file, out_fmt, ref_results)


def test_invalid(make_args, tmp_path):
	"""Test invalid parameter values exit with error code."""

	results_file = tmp_path / ('results.json')

	# No genomes or signatures
	args = make_args(output=results_file)
	assert invoke_cli(args).exit_code != 0

	# Multiple inputs
	args = make_args(output=results_file, positional=True, list_file=True)
	assert invoke_cli(args).exit_code != 0
	args = make_args(output=results_file, positional=True, sig_file=True)
	assert invoke_cli(args).exit_code != 0
	args = make_args(output=results_file, list_file=True, sig_file=True)
	assert invoke_cli(args).exit_code != 0
