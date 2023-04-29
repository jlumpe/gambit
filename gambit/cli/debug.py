from typing import Dict, Any

import click

from .root import cli
from gambit.util.misc import is_importable
from gambit.cli.test import invoke_cli


#: Modules to import in interactive shell.
SHELL_MODULES = dict(
	metric='gambit.metric',
	kmers='gambit.kmers',
)


@cli.group(
	name='debug',
	hidden=True,
)
def debug_group():
	"""Tools for debugging and testing."""
	pass


def make_shell_ns(ctx) -> Dict[str, Any]:
	"""Make the user namespace for the shell command."""
	from importlib import import_module

	ns = dict(
		click_ctx=ctx,
		ctx=ctx.obj,
	)

	# Import modules into namespace
	for alias, name in SHELL_MODULES.items():
		ns[alias] = import_module(name)

	return ns

@debug_group.command()
@click.option(
	'--ipython/--no-ipython', 'use_ipython',
	default=None,
	help='Use IPython instead of built-in Python REPL.',
)
@click.pass_context
def shell(ctx, use_ipython):
	"""Start an interactive shell with application data and modules imported.

	Attempts to launch an IPython interactive interpreter if it is installed,
	otherwise falls back on standard Python REPL.
	"""

	if use_ipython is None:
		if is_importable('IPython'):
			use_ipython = True
		else:
			click.echo('IPython not available, defaulting to built-in Python REPL.', err=True)
			use_ipython = False

	ns = make_shell_ns(ctx)

	if use_ipython:
		from IPython import start_ipython
		start_ipython(argv=[], user_ns=ns)

	else:
		from code import interact
		interact(local=ns)


@debug_group.command()
@click.option(
	'--ipdb/--no-ipdb', 'use_ipdb',
	default=None,
	help='Use IPython debugger instead of pdb.',
)
@click.argument('cmdname', metavar='COMMAND')
@click.argument('args', nargs=-1)
@click.pass_context
def invoke(ctx, cmdname, args, use_ipdb):
	"""Invoke another GAMBIT command within the Python debugger."""

	if use_ipdb is None:
		use_ipdb = is_importable('ipdb')
		if not use_ipdb:
			click.echo('ipdb not installed, defaulting to pdb')

	if use_ipdb:
		import ipdb as db
	else:
		import pdb as db

	try:
		invoke_cli([cmdname, *args])
	except Exception:
		db.post_mortem()
	else:
		click.echo('Command completed successfully')
