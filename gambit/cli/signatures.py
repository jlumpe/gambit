from typing import Optional, List, BinaryIO, TextIO
import sys

import click

from .util import print_table
import gambit.io.json as gjson
from gambit.signatures.hdf5 import HDF5Signatures


def format_none(value):
	return '<none>' if value is None else value


@click.group(name='signatures')
def signatures_group():
	"""Create and inspect GAMBIT signature files."""
	pass


@signatures_group.command()
@click.option(
	'-j', '--json',
	is_flag=True,
	help='Write output in JSON format.',
)
@click.option(
	'-p', '--pretty',
	is_flag=True,
	help='Prettify JSON output.',
)
@click.option(
	'-i', '--ids',
	is_flag=True,
	help='Write IDs of signatures in file, one per line.',
)
@click.argument(
	'file',
	type=click.Path(exists=True, dir_okay=False),
)
def info(file: str, json: bool, pretty: bool, ids: bool):
	"""Inspect GAMBIT signature files."""
	sigs = HDF5Signatures.open(file)

	if ids:
		if json:
			raise click.ClickException('The -i/--ids and -j/--json options are mutually exclusive.')

		for id in sigs.ids:
			click.echo(id)

	elif json:
		data = dict(
			count=len(sigs),
			kmerspec=sigs.kmerspec,
			metadata=sigs.meta,
		)
		kw = dict(indent=' ', sort_keys=True) if pretty else dict()
		gjson.dump(data, sys.stdout, **kw)

	else:
		rows1 = [
			('Count:', len(sigs)),
			('k:', sigs.kmerspec.k),
			('Prefix:', sigs.kmerspec.prefix.decode('ascii')),
			('Data type:', sigs.dtype),
		]
		print_table(rows1, colsep='  ')

		print('Metadata:')

		rows2 = [
			('ID:', format_none(sigs.meta.id)),
			('Name:', format_none(sigs.meta.name)),
			('Version:', format_none(sigs.meta.version )),
			('Description:', format_none(sigs.meta.description)),
			('Genome ID attribute:', format_none(sigs.meta.id_attr)),
			('Has extra:', 'no' if sigs.meta.extra is None else 'yes'),
		]
		print_table(rows2, colsep='  ', left = '  ')
