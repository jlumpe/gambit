"""
Test the 'gambit query' CLI command using the testdb_210818 database.
"""

from copy import copy
from typing import Optional, Iterable
from pathlib import Path

import pytest

from gambit.seq import SequenceFile
from gambit.query import QueryInput, QueryResults
from gambit.util.misc import zip_strict
from gambit.util.io import write_lines, FilePath
from gambit.cli.common import strip_seq_file_ext

from ..testdb import TestDB
from ..results import check_json_results, check_csv_results
from .common import invoke_cli


def make_args(testdb: TestDB, *,
			  positional_files: Optional[Iterable[SequenceFile]] = None,
			  list_file: Optional['FilePath'] = None,
			  sig_file: bool = False,
			  output: Optional['FilePath'] = None,
			  outfmt: Optional[str] = None,
			  strict: bool=False,
			  ) -> list[str]:
	"""Make command line arguments for querying."""

	args: list[str] = [f'--db={testdb.paths.root}', 'query']
	args.append('--strict' if strict else '--no-strict')

	if output is not None:
		args.append(f'--output={output}')

	if outfmt is not None:
		args.append(f'--outfmt={outfmt}')

	if positional_files is not None:
		args.extend(map(str, positional_files))

	if list_file is not None:
		args += ['-l', str(list_file), f'--ldir={testdb.paths.query_genomes_dir}']

	if sig_file:
		args.append(f'--sigfile={testdb.paths.query_signatures}')

	return args


def make_ref_results(testdb: TestDB, inputs: Iterable[QueryInput], strict: bool, nqueries: Optional[int]):
	"""
	Make a copy of the reference query results to compare to, modifying to account for possibly
	different query inputs and # of queries.
	"""
	ref_results = copy(testdb.get_query_results(strict))
	ref_results.items = ref_results.items[:nqueries]

	for item, input in zip_strict(ref_results.items, inputs):
		item.input = input

	return ref_results


def check_results(results_file: Path, out_fmt: str, ref_results: QueryResults):
	"""Check results output matches reference QueryResults object."""
	if out_fmt == 'json':
		with open(results_file) as fh:
			check_json_results(fh, ref_results, strict=False)

	elif out_fmt == 'csv':
		with open(results_file) as fh:
			check_csv_results(fh, ref_results, strict=False)

	elif out_fmt == 'archive':
		assert results_file.is_file()  # TODO

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
		(20,   False, 'archive', False, False),
	],
)
def test_full_query(testdb: TestDB,
					nqueries: Optional[int],
					use_list_file: bool,
					out_fmt: str,
					strict: bool,
					gzipped: bool,
					tmp_path: Path,
					):
	"""Run a full query using the command line interface."""

	query_files = testdb.get_query_files(gzipped)[:nqueries]
	inputs = [
		QueryInput(strip_seq_file_ext(file.path.name), file)
		for file in query_files
	]
	ref_results: QueryResults = make_ref_results(testdb, inputs, strict, nqueries)

	results_file = tmp_path / ('results.' + out_fmt)

	if use_list_file:
		list_file = tmp_path / 'genomes.txt'
		write_lines(query_files, list_file)
		input_kw = dict(list_file=list_file)
	else:
		input_kw = dict(positional_files=query_files)

	args = make_args(
		testdb,
		output=results_file,
		outfmt=out_fmt,
		strict=strict,
		**input_kw,
	)

	invoke_cli(args)
	check_results(results_file, out_fmt, ref_results)


# Not really necessary to check all combinations of parameters.
@pytest.mark.parametrize('out_fmt', ['json'])
@pytest.mark.parametrize('strict', [False])
def test_sigfile(testdb: TestDB, out_fmt: str, strict: bool, tmp_path: Path):
	"""Test using signature file instead of parsing genome files."""

	inputs = list(map(QueryInput, testdb.query_signatures.ids))
	ref_results = make_ref_results(testdb, inputs, strict, None)

	results_file = tmp_path / ('results.' + out_fmt)

	args = make_args(
		testdb,
		sig_file=True,
		output=results_file,
		outfmt=out_fmt,
		strict=False,
	)

	invoke_cli(args)
	check_results(results_file, out_fmt, ref_results)


def test_invalid(testdb: TestDB, tmp_path: Path):
	"""Test invalid parameter values exit with error code."""

	query_files = testdb.get_query_files()
	list_file = tmp_path / 'list.json'
	write_lines(query_files, list_file)
	results_file = tmp_path / ('results.json')

	# No genomes or signatures
	args = make_args(testdb, output=results_file)
	result = invoke_cli(args, success=False)
	assert result.stderr.strip() == 'Error: One of GENOMES, -l, or -s/--sigfile is required'

	# Multiple inputs
	multi_msg = 'Error: GENOMES, -l, and -s/--sigfile are mutually exclusive'

	args = make_args(testdb, output=results_file, positional_files=query_files, list_file=list_file)
	assert invoke_cli(args, success=False).stderr.strip() == multi_msg

	args = make_args(testdb, output=results_file, positional_files=query_files, sig_file=True)
	assert invoke_cli(args, success=False).stderr.strip() == multi_msg

	args = make_args(testdb, output=results_file, list_file=list_file, sig_file=True)
	assert invoke_cli(args, success=False).stderr.strip() == multi_msg
