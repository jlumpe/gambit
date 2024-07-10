"""
Test the 'gambit query' CLI command using the testdb_210818 database.
"""

import os
from copy import copy

import pytest

from gambit.cli.test import invoke_cli
from gambit.results.test import check_json_results, check_csv_results
from gambit.seq import SequenceFile
from gambit.query import QueryInput
from gambit.util.misc import zip_strict
from gambit.util.io import write_lines
from gambit.cli.common import strip_seq_file_ext


@pytest.fixture(params=[None])
def nqueries(request):
	"""Number of testdb query files to use, None means use all of them.

	Can be changed via indirect parameterization in specific tests.
	Note than with slice notation, `[:None]` is the same as `[:]`.
	"""
	return request.param


@pytest.fixture()
def query_files(testdb, nqueries):
	"""Paths to query files."""
	return [SequenceFile(f.path, f.format, f.compression) for f in testdb.get_query_files()[:nqueries]]


@pytest.fixture()
def cd_query_genomes(testdb):
	"""Change working directory to query genomes directory."""
	old_wd = os.getcwd()
	try:
		os.chdir(testdb.paths.query_genomes_dir)
		yield
	finally:
		os.chdir(old_wd)


@pytest.fixture(name='make_args')
def make_args_factory(testdb, query_files, tmp_path):

	def make_args(positional=False, list_file=False, sig_file=False, output=None, outfmt=None, strict=False):
		"""Make command line arguments for query file."""

		args = [f'--db={testdb.paths.root}', 'query']
		args.append('--strict' if strict else '--no-strict')

		if output is not None:
			args.append(f'--output={output}')

		if outfmt is not None:
			args.append(f'--outfmt={outfmt}')

		if positional:
			args.extend(query_files)

		if list_file:
			list_file = tmp_path / 'genomes.txt'
			write_lines(query_files, list_file)
			args += ['-l', str(list_file), f'--ldir={testdb.paths.query_genomes_dir}']

		if sig_file:
			args.append(f'--sigfile={testdb.paths.query_signatures}')

		return list(map(str, args))

	return make_args

@pytest.fixture(name='make_ref_results')
def make_ref_results_factory(testdb, nqueries, query_files):
	"""
	Make a copy of the reference query results to compare to, modifying to account for possibly
	different query inputs and # of queries.
	"""
	def make_ref_results(strict, inputs):
		ref_results = copy(testdb.get_query_results(strict))
		ref_results.items = ref_results.items[:nqueries]

		for item, input in zip_strict(ref_results.items, inputs):
			item.input = input

		return ref_results

	return make_ref_results


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
	['nqueries', 'use_list_file', 'out_fmt', 'strict', 'gzipped'],
	[
		(None, False, 'json', False, False),
		(20,   False, 'csv',  False, False),
		(None, False, 'json', True,  False),
		(20,   False, 'csv',  True,  False),
		(None, False, 'json', False, True),
		(20,   True,  'json', False, False),
	],
	indirect=['nqueries'],
)
def test_full_query(make_args, make_ref_results, use_list_file, out_fmt, strict, gzipped, query_files, tmp_path):
	"""Run a full query using the command line interface."""

	inputs = [
		QueryInput(strip_seq_file_ext(file.path.name), file)
		for file in query_files
	]
	ref_results = make_ref_results(strict, inputs)

	results_file = tmp_path / ('results.' + out_fmt)

	args = make_args(
		positional=not use_list_file,
		list_file=use_list_file,
		output=results_file,
		outfmt=out_fmt,
		strict=strict,
	)

	invoke_cli(args)
	check_results(results_file, out_fmt, ref_results)


# Not really necessary to check all combinations of parameters.
@pytest.mark.parametrize('out_fmt', ['json'])
@pytest.mark.parametrize('strict', [False])
def test_sigfile(make_args, make_ref_results, testdb, out_fmt, strict, tmp_path):
	"""Test using signature file instead of parsing genome files."""

	inputs = list(map(QueryInput, testdb.query_signatures.ids))
	ref_results = make_ref_results(strict, inputs)

	results_file = tmp_path / ('results.' + out_fmt)

	args = make_args(
		sig_file=True,
		output=results_file,
		outfmt=out_fmt,
		strict=False,
	)

	invoke_cli(args)
	check_results(results_file, out_fmt, ref_results)


def test_invalid(make_args, tmp_path):
	"""Test invalid parameter values exit with error code."""

	results_file = tmp_path / ('results.json')

	# No genomes or signatures
	args = make_args(output=results_file)
	invoke_cli(args, success=False)

	# Multiple inputs
	args = make_args(output=results_file, positional=True, list_file=True)
	assert invoke_cli(args, success=False)
	args = make_args(output=results_file, positional=True, sig_file=True)
	assert invoke_cli(args, success=False)
	args = make_args(output=results_file, list_file=True, sig_file=True)
	assert invoke_cli(args, success=False)
