"""Test code in gambit.cli.common."""

from pathlib import Path

import pytest
import click
import numpy as np

from gambit.cli import cli, common
from gambit.cli.test import default_runner, allow_no_args
from gambit.db import ReferenceDatabase
from gambit.seq import SequenceFile
from gambit.util.misc import zip_strict
from gambit.util.io import write_lines


class TestCLIContext:

	def make_context(self, args, env=None):
		runner = default_runner()

		with runner.isolation(env=env), allow_no_args(cli):
			ctx = cli.make_context('gambit', args)

		return common.CLIContext(ctx)

	def test_no_db(self):
		"""Test with no database specified."""
		ctx = self.make_context([])

		assert ctx.db_path is None
		assert not ctx.has_database
		assert not ctx.has_genomes
		assert not ctx.has_signatures
		assert ctx.engine is None
		assert ctx.Session is None
		assert ctx.signatures is None

		with pytest.raises(click.ClickException):
			ctx.require_database()
		with pytest.raises(click.ClickException):
			ctx.require_genomes()
		with pytest.raises(click.ClickException):
			ctx.require_signatures()

	@pytest.mark.parametrize('method', ['option', 'envvar'])
	def test_with_db(self, method, testdb):
		"""Test with database given through the --db argument or environment variable."""
		dbpath = testdb.paths.root

		if method == 'option':
			args = ['--db', str(dbpath)]
			env = dict()
		elif method == 'envvar':
			args = []
			env = dict(GAMBIT_DB_PATH=str(dbpath))
		else:
			assert 0

		ctx = self.make_context(args, env=env)

		assert ctx.db_path == dbpath
		assert ctx.has_database
		assert ctx.has_genomes
		assert ctx.has_signatures
		assert ctx.engine is not None
		assert ctx.Session is not None
		assert ctx.signatures is not None

		ctx.require_database()
		ctx.require_genomes()
		ctx.require_signatures()

		db = ReferenceDatabase.load_from_dir(dbpath)
		ctx_db = ctx.get_database()
		assert db.genomeset.key == ctx_db.genomeset.key
		assert all(g1.key == g2.key for g1, g2 in zip(db.genomes, ctx_db.genomes))
		assert np.array_equal(db.sig_indices, ctx_db.sig_indices)
		assert db.signatures.meta.id == ctx_db.signatures.meta.id
		assert np.array_equal(db.signatures.ids, ctx_db.signatures.ids)


def test_strip_seq_file_ext():
	"""Test the strip_seq_file_ext function."""
	for stem in ['foo', 'GCF_000000000.1']:
		for ext in ['.fasta', '.fna', '.fa']:
			assert common.strip_seq_file_ext(stem) == stem
			assert common.strip_seq_file_ext(stem + ext) == stem
			assert common.strip_seq_file_ext(stem + ext + '.gz') == stem


@pytest.mark.parametrize(
	['strip_dir', 'strip_ext'],
	[(True, True), (True, False), (False, False)],
)
class TestGetSequenceFiles:
	"""Test the get_sequence_files() function."""

	def check_ids(self, ids, paths, strip_dir, strip_ext):
		for id_, path in zip_strict(ids, paths):
			if strip_dir:
				expected = Path(path).name
				if strip_ext:
					expected = expected.split('.')[0]
			else:
				expected = str(path)

			assert id_ == expected

	def check_files(self, files, paths):
		for file, path in zip_strict(files, paths):
			assert isinstance(file, SequenceFile)
			assert file.path == Path(path)
			assert file.format == 'fasta'
			assert file.compression == 'auto'

	def test_explicit(self, strip_dir, strip_ext):
		"""Test given explicit paths from CLI argument."""
		paths = [f'path/to/{i + 1}.fasta' for i in range(10)]
		ids, files = common.get_sequence_files(paths, None, None, strip_dir=strip_dir, strip_ext=strip_ext)
		self.check_ids(ids, paths, strip_dir, strip_ext)
		self.check_files(files, paths)

	@pytest.mark.parametrize('wd,absolute', [
		('.', False),                # Relative to current directory
		('path/to/genomes', False),  # Relative to other directory
		('foo/baz', True),           # Absolute paths in file, ignore wd
	])
	def test_listfile(self, wd, absolute, tmpdir, strip_dir, strip_ext):
		"""Test reading file paths from list file."""
		wd = Path(wd)
		list_paths = [f'{i + 1}.fasta' for i in range(10)]
		if absolute:
			wd2 = Path('/home/jlumpe')
			list_paths = [wd2 / p for p in list_paths]

		listfile = tmpdir / 'files.txt'
		write_lines(list_paths, listfile)
		true_paths = [wd / p for p in list_paths]

		ids, files = common.get_sequence_files([], listfile, wd, strip_dir=strip_dir, strip_ext=strip_ext)
		self.check_ids(ids, list_paths, strip_dir, strip_ext)
		self.check_files(files, true_paths)


def test_params_by_name():
	from gambit.cli.query import query_cmd as cmd

	params = common.params_by_name(cmd)
	assert set(params.values()) == set(cmd.params)
	for name, param in params.items():
		assert param.name == name

	params2 = cmd.params[:3]
	assert common.params_by_name(cmd, [param.name for param in params2]) == params2
