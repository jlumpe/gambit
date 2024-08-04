"""Tests for the "dist" command."""

import json
from typing import Optional, Iterable
from pathlib import Path

import pytest
import numpy as np

from gambit.kmers import KmerSpec
from gambit.metric import jaccarddist_matrix
from gambit.sigs import SignatureList, dump_signatures
from gambit.cli.test import invoke_cli
from gambit.util.io import write_lines
from gambit.cluster import load_dmat_csv
import gambit.util.json as gjson
from gambit.kmers import DEFAULT_KMERSPEC
from gambit.seq import SequenceFile
from gambit.cli.common import strip_seq_file_ext

from ..testdb import TestDB


def get_query_files(testdb: TestDB, n: Optional[int] = None, gz: bool = False) -> list[SequenceFile]:
	return testdb.get_query_files(gz)[:n]


def get_ref_files(testdb: TestDB, n: Optional[int] = None, gz: bool = False) -> list[SequenceFile]:
	return testdb.get_ref_files(gz)[:n]


def make_args(testdb: TestDB,
			  outfile: Path,
			  *,
			  q_opt: Optional[list[SequenceFile]] = None,  # Query files with -q option
              q_list: Optional[Path] = None,               # Query list file
              q_sigs: bool = False,                        # Use query signature file
              r_opt: Optional[list[SequenceFile]] = None,  # Ref files with -r option
              r_list: Optional[Path] = None,               # Ref list file
              r_sigs: bool = False,                        # Use refs signature file
              r_db: bool = False,                          # Use db for refs
              with_db: bool = False,                       # Pass db at root level
              kmerspec: Optional[KmerSpec] = None,         # Pass -k and -p options
              extra: Iterable[str] = (),                   # Additional args
              ) -> list[str]:

	args: list[str] = ['dist', '-o', str(outfile), *extra]

	if with_db:
		args.insert(0, f'--db={testdb.paths.root}')

	# Queries
	if q_opt is not None:
		for file in q_opt:
			args.extend(['-q', str(file)])
	if q_list is not None:
		args.extend(['--ql', str(q_list)])
		args.extend(['--qdir', str(testdb.paths.query_genomes_dir)])
	if q_sigs:
		args.extend(['--qs', str(testdb.paths.query_signatures)])

	# References
	if r_opt is not None:
		for file in r_opt:
			args.extend(['-r', str(file)])
	if r_list is not None:
		args.extend(['--rl', str(r_list)])
		args.extend(['--rdir', str(testdb.paths.ref_genomes_dir)])
	if r_sigs:
		args.extend(['--rs', str(testdb.paths.ref_signatures)])
	if r_db:
		args.append('--use-db')

	if kmerspec is not None:
		args += [
			'-k', str(kmerspec.k),
			'--prefix', kmerspec.prefix_str,
		]

	return args


def check_output(outfile: Path, expected_matrix: np.ndarray, nqueries: Optional[int], nrefs: Optional[int]):
	dmat, row_ids, col_ids = load_dmat_csv(outfile)
	assert np.allclose(dmat, expected_matrix[:nqueries, :nrefs], atol=1e-4)
	# TODO: check row/col IDs


@pytest.fixture(scope='session')
def expected_matrix(testdb: TestDB):
	return jaccarddist_matrix(testdb.query_signatures, testdb.ref_signatures)


@pytest.fixture(scope='session')
def expected_matrix_square(testdb: TestDB):
	return jaccarddist_matrix(testdb.query_signatures, testdb.query_signatures)


@pytest.mark.parametrize(
	'q_type,r_type,queries_gz,refs_gz',
	[
		('sigs', 'sigs', False, False),
		('list', 'sigs', False, False),
		('sigs', 'list', False, False),
		('list', 'list', False, False),
		('opt',  'sigs', False, False),
		('sigs', 'opt',  False, False),
		('sigs', 'db',   False, False),
		('list', 'sigs', True,  False),
		('sigs', 'list', False, True),
	],
)
def test_basic(testdb: TestDB,
			   q_type: str,                  # Query input format
			   r_type: str,                  # Referencer input format
			   queries_gz: bool,             # Use gzipped query files
			   refs_gz: bool,                # Use gzipped reference files
			   expected_matrix: np.ndarray,
			   tmp_path: Path,
			   ):
	"""Test test basic usage, with query/ref sequences/signatures from different sources."""

	# Use only 10 query/reference files if passing by CLI option or by list file
	nqueries = 10 if q_type in ('opt', 'list') else None
	nrefs = 10 if r_type in ('opt', 'list') else None

	outfile = tmp_path / 'out.csv'
	query_files = get_query_files(testdb, nqueries, queries_gz)
	ref_files = get_ref_files(testdb, nrefs, refs_gz)

	# Query sequence specification
	if q_type == 'opt':
		query_kw = dict(q_opt=query_files)
	elif q_type == 'list':
		q_list = tmp_path / 'queries.txt'
		write_lines(query_files, q_list)
		query_kw = dict(q_list=q_list)
	elif q_type == 'sigs':
		query_kw = dict(q_sigs=True)
	else:
		assert False

	# Reference sequence specification
	if r_type == 'opt':
		ref_kw = dict(r_opt=ref_files)
	elif r_type == 'list':
		r_list = tmp_path / 'refs.txt'
		write_lines(ref_files, r_list)
		ref_kw = dict(r_list=r_list)
	elif r_type == 'sigs':
		ref_kw = dict(r_sigs=True)
	elif r_type == 'db':
		ref_kw = dict(r_db=True)
	else:
		assert False

	using_sigfile = q_type == 'sigs' or r_type == 'sigs'

	args = make_args(
		testdb,
		outfile,
		**query_kw,
		**ref_kw,
		kmerspec=None if using_sigfile else testdb.kmerspec,
		with_db=r_type == 'db',
	)
	invoke_cli(args)
	check_output(outfile, expected_matrix, nqueries, nrefs)


def test_default_kspec(testdb: TestDB, tmp_path: Path):
	"""Test that the default KmerSpec is used when not otherwise specified."""

	outfile = tmp_path / 'out.csv'
	q_list = tmp_path / 'queries.txt'
	q_list.touch()
	r_list = tmp_path / 'refs.txt'
	r_list.touch()

	args = make_args(testdb, outfile, q_list=q_list, r_list=r_list, extra=('--dump-params',))

	result = invoke_cli(args)
	params = json.loads(result.stdout)
	assert params['kmerspec'] == gjson.to_json(DEFAULT_KMERSPEC)


def test_kspec_err(testdb: TestDB, tmp_path: Path):
	"""Test selection of k-mer params and errors on inconsistencies."""

	outfile = tmp_path / 'out.csv'

	query_files = get_query_files(testdb, 10)
	query_lf = tmp_path / 'queries.txt'
	write_lines(query_files, query_lf)

	ref_files = get_ref_files(testdb, 10)
	ref_lf = tmp_path / 'refs.txt'
	write_lines(ref_files, ref_lf)

	# Alternate kmerspec
	kspec1 = testdb.kmerspec
	kspec2 = KmerSpec(5, 'AC')
	assert kspec2 != kspec1

	# Create signatures file for alt kspec
	alt_sigfile = tmp_path / 'alt_sigs.gs'
	alt_sigs = SignatureList([], kspec2)
	dump_signatures(alt_sigfile, alt_sigs)

	# Kspec from args inconsistent with query or reference signatures
	msg = (
		f'Error: K-mer search parameters {{}} ({kspec2.k}/{kspec2.prefix_str}) '
		f'do not match those of {{}} ({kspec1.k}/{kspec1.prefix_str}).'
	)

	args = make_args(testdb, outfile, q_sigs=True, r_list=ref_lf, kmerspec=kspec2)
	result = invoke_cli(args, success=False)
	assert result.stderr.strip() == msg.format('from command line options', 'query signatures')

	args = make_args(testdb, outfile, q_list=query_lf, r_sigs=True, kmerspec=kspec2)
	result = invoke_cli(args, success=False)
	assert result.stderr.strip() == msg.format('from command line options', 'reference signatures')

	args = make_args(testdb, outfile, q_list=query_lf, r_db=True, kmerspec=kspec2, with_db=True)
	result = invoke_cli(args, success=False)
	assert result.stderr.strip() == msg.format('from command line options', 'reference signatures')

	# Ref and query signatures have differing kspec
	args = make_args(testdb, outfile, r_sigs=True, extra=('--qs', str(alt_sigfile)))
	result = invoke_cli(args, success=False)
	assert result.stderr.strip() == msg.format('of query signatures', 'reference signatures')

	args = make_args(testdb, outfile, r_db=True, with_db=True, extra=('--qs', str(alt_sigfile)))
	result = invoke_cli(args, success=False)
	assert result.stderr.strip() == msg.format('of query signatures', 'reference signatures')


@pytest.mark.parametrize(
	'q_type,queries_gz',
	[
		('sigs', False),
		('list', False),
		('opt',  False),
		('list', True),
	],
)
def test_square(testdb: TestDB,
				q_type: str,
				queries_gz: bool,
				expected_matrix_square: np.ndarray,
				tmp_path: Path,
				):
	"""Test --square option."""

	outfile = tmp_path / 'out.csv'
	nqueries = 10 if q_type in ('opts', 'list') else None
	query_files = get_query_files(testdb, nqueries, queries_gz)

	# Query sequence specification
	if q_type == 'opt':
		query_kw = dict(q_opt=query_files)
	elif q_type == 'list':
		q_list = tmp_path / 'queries.txt'
		write_lines(query_files, q_list)
		query_kw = dict(q_list=q_list)
	elif q_type == 'sigs':
		query_kw = dict(q_sigs=True)
	else:
		assert False

	args = make_args(
		testdb,
		outfile,
		**query_kw,
		kmerspec=None if q_type == 'sigs' else testdb.kmerspec,
		extra=['--square'],
	)
	invoke_cli(args)

	check_output(outfile, expected_matrix_square, nqueries, nqueries)
