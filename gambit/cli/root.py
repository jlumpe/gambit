"""Define the root CLI command group."""

import click

from gambit import __version__ as GAMBIT_VERSION
from .common import CLIContext, dirpath
from .query import query_cmd
from .signatures import signatures_group
from .debug import debug_group


# Top-level cli group
@click.group()
@click.option(
	'-d', '--db', 'db_path',
	type=dirpath(exists=True),
	envvar='GAMBIT_DB_PATH',
	help='Directory containing GAMBIT database files.',
)
@click.version_option(GAMBIT_VERSION, prog_name='gambit')
@click.pass_context
def cli(ctx: click.Context, **kw):
	"""Tool for rapid taxonomic identification of microbial pathogens from genomic data."""
	ctx.obj = CLIContext(ctx)


# Add sub-commands
cli.add_command(query_cmd)
cli.add_command(signatures_group)
cli.add_command(debug_group)
