"""Tests for the "dist" command."""

import json

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


@pytest.fixture()
def outfile(tmp_path):
	return tmp_path / 'out.csv'

@pytest.fixture(params=[None])
def nqueries(request):
	return request.param

@pytest.fixture(params=[False])
def queries_gz(request):
	return request.param

@pytest.fixture()
def query_files(testdb, nqueries, queries_gz):
	return [f for f in testdb.get_query_files(queries_gz)[:nqueries]]

@pytest.fixture(params=[None])
def nrefs(request):
	return request.param

@pytest.fixture(params=[False])
def refs_gz(request):
	return request.param

@pytest.fixture()
def ref_files(testdb, nrefs, refs_gz):
	return [f for f in testdb.get_ref_files(refs_gz)[:nrefs]]

@pytest.fixture(name='make_args')
def make_args_factory(testdb, query_files, ref_files, outfile, tmp_path):

	def make_args(q_opt=False,       # Pass queries with -q option
	              q_list=False,      # Pass queries with list file
	              q_sigs=False,      # Use query signature file
	              r_opt=False,       # Pass refs with -r option
	              r_list=False,      # Pass refs with list file
	              r_sigs=False,      # Use refs signature file
	              r_db=False,        # Use db for refs
	              with_db=False,     # Pass db at root level
	              with_kspec=False,  # Pass -k and -p options
	              extra=(),          # Additional args
	              ):

		args = ['dist', '-o', outfile, *extra]

		if with_db:
			args.insert(0, f'--db={testdb.paths.root}')

		if q_opt:
			for file in query_files:
				args.extend(['-q', file])
		if q_list:
			qlfile = tmp_path / 'queries.txt'
			write_lines(query_files, qlfile)
			args.extend(['--ql', qlfile])
			args.extend(['--qdir', testdb.paths.query_genomes_dir])
		if q_sigs:
			args.extend(['--qs', testdb.paths.query_signatures])

		if r_opt:
			for file in ref_files:
				args.extend(['-r', file])
		if r_list:
			rlfile = tmp_path / 'refs.txt'
			write_lines(ref_files, rlfile)
			args.extend(['--rl', rlfile])
			args.extend(['--rdir', testdb.paths.ref_genomes_dir])
		if r_sigs:
			args.extend(['--rs', testdb.paths.ref_signatures])
		if r_db:
			args.append('--use-db')

		if with_kspec:
			args += [
				'-k', str(testdb.kmerspec.k),
				f'--prefix={testdb.kmerspec.prefix_str}',
			]

		return args

	return make_args

@pytest.fixture(scope='session')
def expected_matrix(testdb):
	return jaccarddist_matrix(testdb.query_signatures, testdb.ref_signatures)

@pytest.fixture(scope='session')
def expected_matrix_square(testdb):
	return jaccarddist_matrix(testdb.query_signatures, testdb.query_signatures)

@pytest.fixture(name='check_output')
def check_output_factory(outfile, expected_matrix, nqueries, nrefs):
	def check_output():
		dmat, row_ids, col_ids = load_dmat_csv(outfile)
		assert np.allclose(dmat, expected_matrix[:nqueries, :nrefs], atol=1e-4)
		# TODO check row/column IDs

	return check_output


@pytest.mark.parametrize(
	'q_type,r_type,nqueries,nrefs,queries_gz,refs_gz',
	[
		('sigs', 'sigs', None, None, False, False),
		('list', 'sigs', 10,   None, False, False),
		('sigs', 'list', None, 10  , False, False),
		('list', 'list', 10,   10  , False, False),
		('opt',  'sigs', 10,   None, False, False),
		('sigs', 'opt',  None, 10  , False, False),
		('sigs', 'db',   None, None, False, False),
		('list', 'sigs', 10,   None, True,  False),
		('sigs', 'list', None, 10  , False, True),
	],
	indirect=['nqueries', 'nrefs', 'queries_gz', 'refs_gz'],
)
def test_basic(make_args, check_output, q_type, r_type):
	"""Test test basic usage, with query/ref sequences/signatures from different sources."""

	args = make_args(
		q_opt=q_type == 'opt',
		q_list=q_type == 'list',
		q_sigs=q_type == 'sigs',
		r_opt=r_type == 'opt',
		r_list=r_type == 'list',
		r_sigs=r_type == 'sigs',
		r_db=r_type == 'db',
		with_kspec=True,
		with_db=r_type == 'db',
	)
	invoke_cli(args)
	check_output()

def test_kspec(make_args, testdb, tmp_path):
	"""Test selection of k-mer params and errors on inconsistencies."""

	alt_kspec = KmerSpec(6, 'AC')
	assert alt_kspec != testdb.kmerspec
	alt_kspec_args = ['-k', alt_kspec.k, '-p', alt_kspec.prefix_str]

	alt_sigfile = tmp_path / 'alt_sigs.h5'
	alt_sigs = SignatureList([], alt_kspec)
	dump_signatures(alt_sigfile, alt_sigs)

	# Default kspec
	args = make_args(q_list=True, r_list=True, extra=('--dump-params',))
	result = invoke_cli(args)
	params = json.loads(result.stdout)
	assert params['kmerspec'] == gjson.to_json(DEFAULT_KMERSPEC)

	# Kspec from args inconsistent with query or reference signatures
	args = make_args(q_sigs=True, r_sigs=True) + alt_kspec_args
	invoke_cli(args, success=False)
	args = make_args(q_list=True, r_sigs=True) + alt_kspec_args
	invoke_cli(args, success=False)
	args = make_args(q_sigs=True, r_list=True) + alt_kspec_args
	invoke_cli(args, success=False)
	args = make_args(q_sigs=True, r_db=True) + alt_kspec_args
	invoke_cli(args, success=False)

	# Ref and query signatures inconsistent
	args = make_args(r_sigs=True) + ['--qs', alt_sigfile]
	invoke_cli(args, success=False)
	args = make_args(r_db=True) + ['--qs', alt_sigfile]
	invoke_cli(args, success=False)

@pytest.mark.parametrize(
	'q_type,nqueries,queries_gz',
	[
		('sigs', None, False),
		('list', 10,   False),
		('opt',  10,   False),
		('list', 10,   True),
	],
	indirect=['nqueries', 'queries_gz'],
)
def test_square(make_args, q_type, outfile, expected_matrix_square, nqueries):
	"""Test --square option."""

	args = make_args(
		q_opt=q_type == 'opt',
		q_list=q_type == 'list',
		q_sigs=q_type == 'sigs',
		with_kspec=True,
		extra=['--square'],
	)
	invoke_cli(args)

	out_dmat, row_ids, col_ids = load_dmat_csv(outfile)
	assert np.allclose(out_dmat, expected_matrix_square[:nqueries, :nqueries], atol=1e-4)
	assert row_ids == col_ids
