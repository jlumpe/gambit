from typing import Optional, TextIO
import sys

import click
import h5py as h5

from .common import CLIContext, seq_file_params, get_seq_files
from .util import print_table
from gambit.kmers import KmerSpec
import gambit.io.json as gjson
from gambit.signatures import SignaturesMeta, SignatureArray
from gambit.signatures.hdf5 import HDF5Signatures
from gambit.search import calc_file_signatures
from gambit.util.progress import ClickProgressMeter


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


@signatures_group.command()
@seq_file_params()
@click.option(
	'-o', '--output',
	required=True,
	type=click.Path(writable=True),
	help='File path to write to.',
)
@click.option(
	'-p', '--prefix',
	help='K-mer prefix.',
)
@click.option(
	'-k',
	type=int,
	help='Number of nucleotides to recognize AFTER prefix',
)
@click.option(
	'-m', '--meta-json', 'meta_json',
	type=click.File('r'),
	help='JSON file containing metadata to attach to file.',
)
@click.option(
	'-i', '--ids',
	type=click.File('r'),
	help='File containing genome IDs (one per line).',
)
@click.pass_obj
def create(ctxobj: CLIContext,
           output: str,
           prefix: Optional[str],
           k: Optional[int],
           meta_json: Optional[TextIO],
           ids: Optional[TextIO],
           **kw,
           ):
	"""Create k-mer signatures from genome sequences."""

	seqfiles = get_seq_files(kw)

	if prefix is not None or k is not None:
		# KmerSpec from options
		if not (prefix is not None and k is not None):
			raise click.ClickException('Must specify values for both -k and --prefix arguments.')

		kspec = KmerSpec(k, prefix)

	elif ctxobj.has_signatures:
		# KmerSpec from current reference database
		kspec = ctxobj.signatures.kmerspec

	else:
		raise click.ClickException('Must give values for the -k and --prefix options or specify a reference database.')

	if meta_json is not None:
		meta = gjson.load(meta_json, SignaturesMeta)
	else:
		meta = None

	if ids is not None:
		ids = [line.strip() for line in ids.readlines()]
		if len(ids) != len(seqfiles):
			raise click.ClickException(f'Number of IDs ({len(ids)}) does not match number of genomes ({len(seqfiles)}).')

	else:
		ids = [f.path.name for f in seqfiles]

	sigs = calc_file_signatures(kspec, seqfiles, progress=ClickProgressMeter)
	sigs = SignatureArray(sigs, dtype=kspec.index_dtype)

	with h5.File(output, 'w') as f:
		HDF5Signatures.create(f, kspec, sigs, ids, meta)
