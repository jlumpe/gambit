from typing import Optional, TextIO, List
import sys

import click

from .common import CLIContext, genome_files_arg, print_table, filepath, dirpath, kspec_params, \
	kspec_from_params, read_genomes_list_file
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
	type=filepath(exists=True),
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
			('Has extra:', 'yes' if sigs.meta.extra else 'no'),
		]
		print_table(rows2, colsep='  ', left='  ')


@signatures_group.command()
@genome_files_arg()
@click.option('-l', type=click.File('r'), help='File containing paths to genomes.')
@click.option('--ldir', type=dirpath(), default='.', help='Parent directory of paths in -l.')
@kspec_params
@click.option(
	'-o', '--output',
	required=True,
	type=filepath(writable=True),
	help='File path to write to.',
)
@click.option(
	'-m', '--meta-json', 'meta_file',
	type=click.File('r'),
	help='JSON file containing metadata to attach.',
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
           l: Optional[TextIO],
           ldir: Optional[str],
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

	# Get sequence files
	if files:
		if l is not None:
			raise click.ClickException('-l and GENOMES are mutually exclusive.')
	elif l is not None:
		with l:
			files = read_genomes_list_file(l, ldir)
	else:
		raise click.ClickException('Must give value(s) for either -l or GENOMES.')

	seqfiles = SequenceFile.from_paths(files, 'fasta', 'auto')

	# Get kmerspec
	kspec = kspec_from_params(k, prefix)

	if db_params:
		if kspec is None:
			ctxobj.require_signatures()
			kspec = ctxobj.signatures.kmerspec
		else:
			raise click.ClickException('The -k/--prefix and --db-params options are mutually exclusive.')

	elif kspec is None:
		raise click.ClickException('Must give values for -k/--prefix or specify --db-params.')

	# Metadata / IDs
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

	# Calculate and save
	sigs = calc_file_signatures(kspec, seqfiles, progress=ClickProgressMeter)
	sigs = AnnotatedSignatures(sigs, ids, meta)
	dump_signatures(output, sigs)
