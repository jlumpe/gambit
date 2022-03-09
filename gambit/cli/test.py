"""Tools for testing CLI."""

from typing import Optional, ContextManager, Sequence
from contextlib import contextmanager

import click
from click.testing import CliRunner, Result

from .root import cli


DEFAULT_ENV = dict(
	GAMBIT_DB_PATH=None,  # Ensure this is unset by default in tests.
)


def pop_kwargs(d, keys):
	out = dict()
	for k in keys:
		try:
			v = d.pop(k)
		except KeyError:
			continue
		else:
			out[k] = v

	return out


def default_runner(**kw) -> CliRunner:
	"""Get a CliRunner instance with altered default settings."""
	kw.setdefault('mix_stderr', False)
	kw.setdefault('env', DEFAULT_ENV)
	return CliRunner(**kw)

def invoke_cli(args: Sequence, runner: Optional[CliRunner] = None, **kw) -> Result:
	"""Invoke CLI in test context, using different defaults than base Click method."""
	if runner is None:
		runner = default_runner(**pop_kwargs(kw, ['charset', 'echo_stdin', 'mix_stderr']))

	kw.setdefault('catch_exceptions', False)
	args = list(map(str, args))
	return runner.invoke(cli, args, **kw)


@contextmanager
def allow_no_args(command: click.Command) -> ContextManager[click.Command]:
	"""Context manager which patches a command to allow calling with no arguments.

	Group commands will print help and exit if called without a subcommand, this will also happen
	when just calling the :meth:`click.BaseCommand.create_context` method which interferes with
	testing the context. This function sets the :attr:`click.BaseCommand.no_args_is_help` attribute
	to `False` to disable this behavior, and restores the original value on exit.
	"""
	old_naih = command.no_args_is_help

	try:
		command.no_args_is_help = False
		yield command

	finally:
		command.no_args_is_help = old_naih

