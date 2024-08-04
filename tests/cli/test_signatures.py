"""Tests for the "signatures" command group."""

import json
from pathlib import Path

import pytest
import numpy as np

import gambit.util.json as gjson
from gambit.sigs import SignaturesMeta, load_signatures
from gambit.util.io import write_lines
from gambit.cli.common import strip_seq_file_ext
from gambit.kmers import DEFAULT_KMERSPEC

from ..testdb import TestDB
from .common import invoke_cli


class TestInfoCommand:

	@pytest.fixture(params=[False, True])
	def use_db(self, request):
		return request.param

	@pytest.fixture()
	def base_args(self, testdb: TestDB, use_db: bool):
		if use_db:
			return [f'--db={testdb.paths.root}', 'signatures', 'info', '-d']
		else:
			return ['signatures', 'info', str(testdb.paths.ref_signatures)]

	def test_standard(self, base_args: list[str]):
		result = invoke_cli(base_args)

		# TODO: check

	def test_json(self, base_args: list[str], testdb: TestDB):
		args = [*base_args, '--json']
		result = invoke_cli(args)

		data = json.loads(result.stdout)
		assert data['count'] == len(testdb.ref_signatures)
		assert data['kmerspec'] == gjson.to_json(testdb.ref_signatures.kmerspec)
		assert data['metadata'] == gjson.to_json(testdb.ref_signatures.meta)

	def test_ids(self, base_args: list[str], testdb: TestDB):
		args = [*base_args, '-i']
		result = invoke_cli(args)

		assert np.array_equal(result.stdout.splitlines(), testdb.ref_signatures.ids)

	def test_invalid(self, testdb: TestDB):
		args = [
			f'--db={testdb.paths.root}',
			'signatures',
			'info',
			'-d',
			testdb.paths.ref_signatures,
		]
		invoke_cli(args, success=False)


class TestCreateCommand:

	@pytest.fixture(params=[False])
	def infiles(self, request, testdb: TestDB):
		"""Input files. Parameter is whether or not they are gzipped."""
		return [f.path for f in testdb.get_query_files(request.param)]

	@pytest.fixture()
	def outfile(self, tmp_path: Path):
		return tmp_path / 'signatures.gs'

	@pytest.fixture(name='make_args')
	def make_args_factory(self, outfile: Path, testdb: TestDB, infiles: list[Path], tmp_path: Path):

		def make_args(opts=(), root_args=(), with_kspec=True, positional_files=True, list_file=False):
			args = list(root_args)
			args += [
				'signatures', 'create',
				f'--output={outfile}',
			]
			if with_kspec:
				args += [
					'-k', str(testdb.kmerspec.k),
					f'--prefix={testdb.kmerspec.prefix_str}',
				]
			args.extend(opts)

			if positional_files:
				args += infiles
			if list_file:
				list_file = tmp_path / 'input-files.txt'
				write_lines([f.name for f in infiles], list_file)
				args += ['-l', str(list_file), f'--ldir={testdb.paths.query_genomes_dir}']

			return list(map(str, args))

		return make_args

	@pytest.fixture()
	def default_ids(self, infiles: list[Path]):
		return [strip_seq_file_ext(file.name) for file in infiles]

	@pytest.fixture(name='check_output')
	def check_output_factory(self, outfile: Path, testdb: TestDB, infiles: list[Path], default_ids: list[str]):

		def check_output(expected_ids=default_ids):
			out = load_signatures(outfile)
			assert out == testdb.query_signatures  # Checks contents and .kmerspec
			assert np.array_equal(out.ids, expected_ids)

			return out

		return check_output

	@pytest.mark.parametrize('infiles', [False, True], indirect=True)
	def test_basic(self, make_args, check_output, infiles):
		"""Test with basic arguments."""
		args = make_args()
		invoke_cli(args)
		check_output()

	def test_list_file(self, make_args, infiles: list[Path], default_ids: list[str]):
		"""Test getting genome list from file."""

		args = make_args(['--dump-params'], positional_files=False, list_file=True)
		result = invoke_cli(args)
		params = json.loads(result.stdout)
		assert params['files'] == list(map(str, infiles))
		assert params['ids'] == default_ids

	def test_with_metadata(self, testdb: TestDB, make_args, check_output, tmp_path: Path):
		"""Test with ids and metadata JSON added."""
		# Metadata file
		metadata = SignaturesMeta(
			name='Test signatures',
			version='0.0',
			description='Signatures of some testdb query genomes.',
			extra=dict(foo=3),
		)

		meta_file = tmp_path / 'metadata.json'
		with open(meta_file, 'w') as f:
			gjson.dump(metadata, f)

		# IDs file
		ids = [f'seq-{i}' for i in range(len(testdb.query_genomes))]
		id_file = tmp_path / 'ids.txt'
		write_lines(ids, id_file)

		# Run
		args = make_args([
			f'--ids={id_file}',
			f'--meta-json={meta_file}',
		])

		invoke_cli(args)

		out = check_output(ids)
		assert out.meta == metadata

	def test_kspec_from_refdb(self, make_args, testdb: TestDB):
		"""Test with KmerSpec taken from reference database."""
		args = make_args(
			['-d', '--dump-params'],
			[f'--db={testdb.paths.root}'],
			with_kspec=False,
		)
		result = invoke_cli(args)
		params = json.loads(result.stdout)
		assert params['kmerspec'] == gjson.to_json(testdb.kmerspec)

	def test_default_kspec(self, make_args, testdb: TestDB):
		"""Test with default KmerSpec."""
		args = make_args(
			['--dump-params'],
			with_kspec=False,
		)
		result = invoke_cli(args)
		params = json.loads(result.stdout)
		assert params['kmerspec'] == gjson.to_json(DEFAULT_KMERSPEC)

	def test_invalid(self, testdb: TestDB, make_args):
		"""Test with invalid parameter combinations."""

		# No genomes
		args = make_args(positional_files=False, list_file=False)
		invoke_cli(args, success=False)

		# Positional args and list file
		args = make_args(positional_files=True, list_file=True)
		invoke_cli(args, success=False)

		# Only -k/--prefix
		args = make_args(['-k', str(testdb.kmerspec.k)], with_kspec=False)
		invoke_cli(args, success=False)
		args = make_args(['--prefix', testdb.kmerspec.prefix_str], with_kspec=False)
		invoke_cli(args, success=False)

		# Both both plus --db-params
		args = make_args(['-d'])
		invoke_cli(args, success=False)

		# --db-params without specifying db in root command group
		args = make_args(['-d'], with_kspec=False)
		invoke_cli(args, success=False)

	def test_ids_wrong_len(self, testdb: TestDB, make_args, tmp_path: Path):
		"""Test where number of IDs does not match query files."""

		ids = [f'seq-{i}' for i in range(len(testdb.query_genomes) - 1)]
		id_file = tmp_path / 'ids2.txt'
		write_lines(ids, id_file)

		args = make_args(['--ids', str(id_file)])
		invoke_cli(args, success=False)
