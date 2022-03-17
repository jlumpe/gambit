"""Test code in gambit.cli.common."""

from pathlib import Path

import pytest
import click
import numpy as np

from gambit.cli import cli, common
from gambit.cli.test import default_runner, allow_no_args
from gambit.db import ReferenceDatabase
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


class TestGetSequenceFiles:
	"""Test the get_sequence_files() function."""

	def test_explicit(self):
		"""Test given explicit paths from CLI argument."""
		paths = [f'path/to/{i + 1}.fasta' for i in range(10)]
		ids, paths2 = common.get_sequence_files(paths, None, None)
		assert ids == paths
		assert paths2 == list(map(Path, paths))

	@pytest.mark.parametrize('wd,absolute', [
		('.', False),                # Relative to current directory
		('path/to/genomes', False),  # Relative to other directory
		('.', True),                 # Absolute paths in file, ignore wd
	])
	def test_listfile(self, wd, absolute, tmpdir):
		"""Test reading file paths from list file."""
		paths = [f'{i + 1}.fasta' for i in range(10)]
		if absolute:
			wd2 = '/home/jlumpe'
			paths = [f'{wd2}/{name}' for name in paths]

		listfile = tmpdir / 'files.txt'
		write_lines(paths, listfile)

		ids, paths2 = common.get_sequence_files([], listfile, wd)

		for p1, id, p2 in zip_strict(paths, ids, paths2):
			assert id == p1
			if absolute:
				assert p2 == Path(p1)
			else:
				assert p2 == Path(wd) / p1


def test_params_by_name():
	from gambit.cli.query import query_cmd as cmd

	params = common.params_by_name(cmd)
	assert set(params.values()) == set(cmd.params)
	for name, param in params.items():
		assert param.name == name

	params2 = cmd.params[:3]
	assert common.params_by_name(cmd, [param.name for param in params2]) == params2
