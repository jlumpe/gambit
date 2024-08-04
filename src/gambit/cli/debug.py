from typing import Any, Optional

import click

from .root import cli


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


def make_shell_ns(ctx) -> dict[str, Any]:
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
def shell(ctx, use_ipython: Optional[bool]):
	"""Start an interactive shell with application data and modules imported.

	Attempts to launch an IPython interactive interpreter if it is installed,
	otherwise falls back on standard Python REPL.
	"""
	if use_ipython is None:
		try:
			import IPython
			use_ipython = True
		except ImportError:
			click.echo('IPython not available, defaulting to built-in Python REPL.', err=True)
			use_ipython = False

	ns = make_shell_ns(ctx)

	if use_ipython:
		from IPython import start_ipython
		start_ipython(argv=[], user_ns=ns)

	else:
		from code import interact
		interact(local=ns)
