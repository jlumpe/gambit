"""Tools for testing CLI."""

from typing import Optional

from click.testing import CliRunner, Result

from .root import cli


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
	return CliRunner(**kw)

def invoke_cli(args, runner: Optional[CliRunner] = None, **kw) -> Result:
	"""Invoke CLI in test context, using different defaults than base Click method."""
	if runner is None:
		runner = default_runner(**pop_kwargs(kw, ['charset', 'echo_stdin', 'mix_stderr']))

	kw.setdefault('catch_exceptions', False)
	return runner.invoke(cli, args, **kw)

