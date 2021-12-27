"""Test code in gambit.cli.common."""

import pytest
import click
import numpy as np

from gambit.cli import cli
from gambit.cli.common import CLIContext
from gambit.cli.test import default_runner, allow_no_args
from gambit.db import load_db_from_dir


class TestCLIContext:

	def make_context(self, args, env=None):
		runner = default_runner()

		with runner.isolation(env=env), allow_no_args(cli):
			ctx = cli.make_context('gambit', args)

		return CLIContext(ctx)

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
	def test_with_db(self, method, testdb_files):
		"""Test with database given through the --db argument or environment variable."""
		dbpath = testdb_files['root']

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

		db = load_db_from_dir(dbpath)
		ctx_db = ctx.get_database()
		assert db.genomeset.key == ctx_db.genomeset.key
		assert all(g1.key == g2.key for g1, g2 in zip(db.genomes, ctx_db.genomes))
		assert np.array_equal(db.sig_indices, ctx_db.sig_indices)
		assert db.signatures.meta.id == ctx_db.signatures.meta.id
		assert np.array_equal(db.signatures.ids, ctx_db.signatures.ids)
