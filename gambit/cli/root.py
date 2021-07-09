"""Define the root CLI command group."""

import click

from gambit import __version__ as GAMBIT_VERSION
from .context import CLIContext


# Top-level cli group
@click.group()
@click.option(
	'-d', '--db', 'db_path',
	type=click.Path(exists=True, file_okay=False),
	envvar='GAMBIT_DB_PATH',
	help='Directory containing GAMBIT database files.'
)
@click.version_option(GAMBIT_VERSION, prog_name='gambit')
@click.pass_context
def cli(ctx: click.Context, db_path):
	"""Tool for rapid taxonomic identification of microbial pathogens from genomic data."""
	ctx.obj = CLIContext(db_path)


# Add sub-commands
from .query import query
cli.add_command(query)
