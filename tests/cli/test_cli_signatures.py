"""Tests for the "signatures" command group."""

import json

import pytest
import numpy as np

from gambit.cli.test import invoke_cli
import gambit.io.json as gjson
from gambit.signatures import SignaturesMeta
from gambit.signatures.hdf5 import HDF5Signatures


class TestInfoCommand:

	@pytest.fixture()
	def base_args(self, testdb_files):
		return ['signatures', 'info', str(testdb_files['ref_signatures'])]

	def test_standard(self, base_args, testdb_signatures):
		result = invoke_cli(base_args)
		assert result.exit_code == 0

		# TODO: check

	@pytest.mark.parametrize('use_ref_sigs', [False, True])
	def test_json(self, use_ref_sigs, testdb_files, testdb_signatures):
		if use_ref_sigs:
			args = [f'--db={testdb_files["root"]}', 'signatures', 'info', '--json']
		else:
			args = ['signatures', 'info', str(testdb_files['ref_signatures']), '--json']

		result = invoke_cli(args)
		assert result.exit_code == 0

		data = json.loads(result.stdout)
		assert data['count'] == len(testdb_signatures)
		assert data['kmerspec'] == gjson.to_json(testdb_signatures.kmerspec)
		assert data['metadata'] == gjson.to_json(testdb_signatures.meta)

	def test_ids(self, base_args, testdb_signatures):
		result = invoke_cli([*base_args, '-i'])
		assert result.exit_code == 0

		assert np.array_equal(result.stdout.splitlines(), testdb_signatures.ids)


class TestCreateCommand:
	NGENOMES = 10
	IDS = [f'seq-{i}' for i in range(NGENOMES)]

	@pytest.fixture()
	def kspec(self, testdb_query_signatures):
		return testdb_query_signatures.kmerspec

	@pytest.fixture()
	def seq_files(self, testdb_queries):
		return [q['file'] for q in testdb_queries[:self.NGENOMES]]

	@pytest.fixture()
	def outfile(self, tmp_path):
		return tmp_path / 'signatures.h5'

	@pytest.fixture()
	def id_file(self, tmp_path):
		p = tmp_path / 'ids.txt'
		with open(p, 'w') as f:
			f.writelines(id + '\n' for id in self.IDS)
		return p

	@pytest.fixture(name='make_args')
	def make_args_factory(self, outfile, seq_files):

		def make_args(opts, root_args=None):
			args = [] if root_args is None else root_args
			args += ['signatures', 'create']
			args += [
				f'--output={outfile}',
				'--seqfmt=fasta',
			]
			args += opts
			args += [str(f.path) for f in seq_files]
			return args

		return make_args

	@pytest.fixture(name='check_output')
	def check_output_factory(self, outfile, seq_files, testdb_query_signatures):

		def check_output(default_ids=True):
			out = HDF5Signatures.open(outfile)

			assert out.kmerspec == testdb_query_signatures.kmerspec
			assert out[:] == testdb_query_signatures[:self.NGENOMES]

			if default_ids:
				expected = [f.path.name for f in seq_files]
				assert np.array_equal(out.ids, expected)

			return out

		return check_output

	def test_basic(self, kspec, make_args, check_output):
		"""Test with basic arguments."""
		args = make_args([
			'-k', str(kspec.k),
			f'--prefix={kspec.prefix_str}',
		])

		result = invoke_cli(args)
		assert result.exit_code == 0

		check_output()

	def test_with_metadata(self, kspec, make_args, check_output, id_file, tmp_path):
		"""Test with ids and metadata JSON added."""
		metadata = SignaturesMeta(
			name='Test signatures',
			version='0.0',
			description='Signatures of some testdb query genomes.',
			extra=dict(foo=3),
		)

		meta_file = tmp_path / 'metadata.json'
		with open(meta_file, 'w') as f:
			gjson.dump(metadata, f)

		args = make_args([
			'-k', str(kspec.k),
			f'--prefix={kspec.prefix_str}',
			f'--ids={id_file}',
			f'--meta-json={meta_file}',
		])

		result = invoke_cli(args)
		assert result.exit_code == 0

		out = check_output(default_ids=False)
		assert np.array_equal(out.ids, self.IDS)
		assert out.meta == metadata

	def test_kspec_from_refdb(self, make_args, check_output, testdb_files):
		"""Test with KmerSpec taken from reference database."""
		args = make_args([], root_args=[f'--db={testdb_files["root"]}'])

		result = invoke_cli(args)
		assert result.exit_code == 0

		check_output()

	def test_bad_kspec(self, kspec, make_args):
		"""Test with KmerSpec incorrectly specified."""

		# No -k, --prefix, or --db
		args = make_args([])
		result = invoke_cli(args)
		assert result.exit_code != 0

		# Only -k
		args = make_args(['-k', str(kspec.k)])
		result = invoke_cli(args)
		assert result.exit_code != 0

		# Only --prefix
		args = make_args(['--prefix', kspec.prefix_str])
		result = invoke_cli(args)
		assert result.exit_code != 0

	def test_ids_wrong_len(self, kspec, make_args, tmp_path):
		"""Test number of IDs do not match query files."""

		id_file = tmp_path / 'ids2.txt'
		with open(id_file, 'w') as f:
			f.writelines(id + '\n' for id in self.IDS[:-1])

		args = make_args([
			'-k', str(kspec.k),
			'--prefix', kspec.prefix_str,
			'--ids', str(id_file)
		])

		result = invoke_cli(args)
		assert result.exit_code != 0
