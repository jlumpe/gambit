from typing import Optional, TextIO, List
import sys

import click

from .common import CLIContext, genome_files_arg, print_table
from gambit.kmers import KmerSpec
import gambit.util.json as gjson
from gambit.sigs import SignaturesMeta, AnnotatedSignatures, load_signatures, dump_signatures
from gambit.sigs.calc import calc_file_signatures
from gambit.util.progress import ClickProgressMeter
from gambit.seq import SequenceFile


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
@click.option(
	'-d', 'use_db',
	is_flag=True,
	help='Use signatures from reference database.',
)
@click.argument(
	'file',
	type=click.Path(exists=True, dir_okay=False),
	required=False,
)
@click.pass_obj
def info(ctxobj: CLIContext, file: str, json: bool, pretty: bool, ids: bool, use_db: bool):
	"""Inspect GAMBIT signature files."""

	if use_db == (file is not None):
		raise click.ClickException('Must specify exactly one of FILE or -d')
	elif file is not None:
		sigs = load_signatures(file)
	else:
		ctxobj.require_signatures()
		sigs = ctxobj.signatures

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
@genome_files_arg()
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
	'-m', '--meta-json', 'meta_file',
	type=click.File('r'),
	help='JSON file containing metadata to attach to file.',
)
@click.option(
	'-i', '--ids', 'ids_file',
	type=click.File('r'),
	help='File containing genome IDs (one per line).',
)
@click.option(
	'-d', '--db-params',
	is_flag=True,
	help='Use k/prefix from reference database.'
)
# Dump parsed CLI parameters and exit. For testing.
@click.option('--dump-params', is_flag=True, hidden=True)
@click.pass_obj
def create(ctxobj: CLIContext,
           files: List[str],
           output: str,
           prefix: Optional[str],
           k: Optional[int],
           meta_file: Optional[TextIO],
           ids_file: Optional[TextIO],
           db_params: bool,
           dump_params: bool,
           ):
	"""Create k-mer signatures from genome sequences."""

	seqfiles = SequenceFile.from_paths(files, 'fasta', 'auto')

	# Get kmerspec
	if prefix is not None or k is not None:
		# KmerSpec from options
		if prefix is None or k is None:
			raise click.ClickException('Must specify values for both -k and --prefix arguments.')
		if db_params:
			raise click.ClickException('The -k/--prefix and -d/--db-params options are mutually exclusive.')

		kspec = KmerSpec(k, prefix)

	elif db_params:
		# KmerSpec from current reference database
		ctxobj.require_signatures()
		kspec = ctxobj.signatures.kmerspec

	else:
		raise click.ClickException('Must give values for the -k and --prefix options or specify a reference database.')

	# Get metadata
	if meta_file is not None:
		meta = gjson.load(meta_file, SignaturesMeta)
	else:
		meta = None

	if ids_file is not None:
		ids = [line.strip() for line in ids_file.readlines()]
		if len(ids) != len(seqfiles):
			raise click.ClickException(f'Number of IDs ({len(ids)}) does not match number of genomes ({len(seqfiles)}).')

	else:
		ids = [f.path.name for f in seqfiles]

	# Dump parameters
	if dump_params:
		params = dict(
			kmerspec=kspec,
			files=files,
			meta=meta,
			ids=ids,
		)
		gjson.dump(params, sys.stdout)
		return

	sigs = calc_file_signatures(kspec, seqfiles, progress=ClickProgressMeter)
	sigs = AnnotatedSignatures(sigs, ids, meta)

	dump_signatures(output, sigs)
