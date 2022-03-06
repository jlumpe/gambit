"""Tests for the "signatures" command group."""

import json

import pytest
import numpy as np

from gambit.cli.test import invoke_cli
import gambit.util.json as gjson
from gambit.sigs import SignaturesMeta, load_signatures
from gambit.util.io import write_lines


class TestInfoCommand:

	@pytest.fixture(params=[False, True])
	def use_db(self, request):
		return request.param

	@pytest.fixture()
	def base_args(self, testdb_files, use_db):
		if use_db:
			return [f'--db={testdb_files["root"]}', 'signatures', 'info', '-d']
		else:
			return ['signatures', 'info', str(testdb_files['ref_signatures'])]

	def test_standard(self, base_args, testdb_signatures):
		result = invoke_cli(base_args)
		assert result.exit_code == 0

		# TODO: check

	def test_json(self, base_args, testdb_signatures):
		args = [*base_args, '--json']
		result = invoke_cli(args)
		assert result.exit_code == 0

		data = json.loads(result.stdout)
		assert data['count'] == len(testdb_signatures)
		assert data['kmerspec'] == gjson.to_json(testdb_signatures.kmerspec)
		assert data['metadata'] == gjson.to_json(testdb_signatures.meta)

	def test_ids(self, base_args, testdb_files, testdb_signatures):
		args = [*base_args, '-i']
		result = invoke_cli(args)
		assert result.exit_code == 0

		assert np.array_equal(result.stdout.splitlines(), testdb_signatures.ids)

	def test_invalid(self, testdb_files):
		assert invoke_cli(['signatures', 'info']).exit_code != 0

		args = [
			f'--db={testdb_files["root"]}',
			'signatures',
			'info',
			'-d',
			str(testdb_files['ref_signatures']),
		]
		assert invoke_cli(args).exit_code != 0


class TestCreateCommand:
	NGENOMES = 10

	@pytest.fixture()
	def kspec(self, testdb_query_signatures):
		return testdb_query_signatures.kmerspec

	@pytest.fixture()
	def outfile(self, tmp_path):
		return tmp_path / 'signatures.h5'

	@pytest.fixture()
	def seq_files(self, testdb_query_files):
		return testdb_query_files[:self.NGENOMES]

	@pytest.fixture(name='make_args')
	def make_args_factory(self, outfile, seq_files, kspec, tmp_path):

		def make_args(opts=(), root_args=(), with_kspec=True, positional_files=True, list_file=False):
			args = list(root_args)
			args += [
				'signatures', 'create',
				f'--output={outfile}',
			]
			if with_kspec:
				args += [
					'-k', str(kspec.k),
					f'--prefix={kspec.prefix_str}',
				]
			args.extend(opts)

			if positional_files:
				args += [str(f.path) for f in seq_files]
			if list_file:
				list_file = tmp_path / 'input-files.txt'
				write_lines([f.path.name for f in seq_files], list_file)
				args += ['-l', str(list_file), f'--ldir={seq_files[0].path.parent}',]

			return args

		return make_args

	@pytest.fixture(name='check_output')
	def check_output_factory(self, outfile, seq_files, testdb_query_signatures):

		def check_output(ids=None):
			out = load_signatures(outfile)

			assert out.kmerspec == testdb_query_signatures.kmerspec
			assert out == testdb_query_signatures[:self.NGENOMES]

			if ids is None:
				ids = [f.path.name for f in seq_files]
			assert np.array_equal(out.ids, ids)

			return out

		return check_output

	@pytest.mark.parametrize('testdb_queries_gzipped', [False, True], indirect=True)
	def test_basic(self, make_args, check_output):
		"""Test with basic arguments."""
		args = make_args()
		result = invoke_cli(args)
		assert result.exit_code == 0

		check_output()

	def test_list_file(self, make_args, seq_files):
		"""Test getting genome list from file."""

		args = make_args(['--dump-params'], positional_files=False, list_file=True)
		result = invoke_cli(args)
		assert result.exit_code == 0
		params = json.loads(result.stdout)
		assert params['files'] == [str(f.path) for f in seq_files]

	def test_with_metadata(self, make_args, check_output, tmp_path):
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
		ids = [f'seq-{i}' for i in range(self.NGENOMES)]
		id_file = tmp_path / 'ids.txt'
		write_lines(ids, id_file)

		# Run
		args = make_args([
			f'--ids={id_file}',
			f'--meta-json={meta_file}',
		])

		result = invoke_cli(args)
		assert result.exit_code == 0

		out = check_output(ids)
		assert out.meta == metadata

	def test_kspec_from_refdb(self, make_args, testdb_files, testdb_signatures):
		"""Test with KmerSpec taken from reference database."""
		args = make_args(
			['-d', '--dump-params'],
			[f'--db={testdb_files["root"]}'],
			with_kspec=False,
		)
		result = invoke_cli(args)
		assert result.exit_code == 0
		params = json.loads(result.stdout)
		assert params['kmerspec'] == gjson.to_json(testdb_signatures.kmerspec)

		# Without specifying db in root command group
		args = make_args(['-d', '--dump-params'], with_kspec=False)
		result = invoke_cli(args)
		assert result.exit_code != 0

	def test_invalid(self, kspec, make_args):
		"""Test with invalid parameter combinations."""

		# No genomes
		args = make_args(positional_files=False, list_file=False)
		assert invoke_cli(args).exit_code != 0

		# Positional args and list file
		args = make_args(positional_files=True, list_file=True)
		assert invoke_cli(args).exit_code != 0

		# No -k, --prefix, or --db-params
		args = make_args(with_kspec=False)
		assert invoke_cli(args).exit_code != 0

		# Only -k/--prefix
		args = make_args(['-k', str(kspec.k)], with_kspec=False)
		assert invoke_cli(args).exit_code != 0
		args = make_args(['--prefix', kspec.prefix_str], with_kspec=False)
		assert invoke_cli(args).exit_code != 0

		# Both both plus --db-params
		args = make_args(['-d'])
		assert invoke_cli(args).exit_code != 0

	def test_ids_wrong_len(self, kspec, make_args, tmp_path):
		"""Test number of IDs do not match query files."""

		ids = [f'seq-{i}' for i in range(self.NGENOMES - 1)]
		id_file = tmp_path / 'ids2.txt'
		write_lines(ids, id_file)

		args = make_args(['--ids', str(id_file)])
		result = invoke_cli(args)
		assert result.exit_code != 0
