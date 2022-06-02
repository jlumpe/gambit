from typing import Optional, TextIO, List
import sys

import click

from . import common
from .root import cli
import gambit.util.json as gjson
from gambit.sigs import SignaturesMeta, AnnotatedSignatures, load_signatures, dump_signatures
from gambit.sigs.calc import calc_file_signatures
from gambit.seq import SequenceFile


def format_none(value):
	return '<none>' if value is None else value


@cli.group(name='signatures')
def signatures_group():
	"""Create and inspect GAMBIT signature files."""
	pass


@signatures_group.command(no_args_is_help=True)
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
	type=common.filepath(exists=True),
	required=False,
)
@click.pass_context
def info(ctx: click.Context, file: str, json: bool, pretty: bool, ids: bool, use_db: bool):
	"""Inspect GAMBIT signature files."""

	common.check_params_group(ctx, ['file', 'use_db'], True, True)
	common.check_params_group(ctx, ['ids', 'json'], True, False)

	if file is not None:
		sigs = load_signatures(file)
	elif use_db:
		ctxobj = ctx.obj  # type: common.CLIContext
		ctxobj.require_signatures()
		sigs = ctx.obj.signatures
	else:
		assert 0

	if ids:
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
		common.print_table(rows1, colsep='  ')

		print('Metadata:')

		rows2 = [
			('ID:', format_none(sigs.meta.id)),
			('Name:', format_none(sigs.meta.name)),
			('Version:', format_none(sigs.meta.version )),
			('Description:', format_none(sigs.meta.description)),
			('Genome ID attribute:', format_none(sigs.meta.id_attr)),
			('Has extra:', 'yes' if sigs.meta.extra else 'no'),
		]
		common.print_table(rows2, colsep='  ', left='  ')


@signatures_group.command(no_args_is_help=True)
@common.genome_files_arg()
@click.option(
	'-l', 'list_file',
	type=click.File('r'),
	metavar='LISTFILE',
	help='File containing names/paths of genome files.',
)
@click.option('--ldir', type=common.dirpath(), default='.', help='Parent directory of paths in LISTFILE.')
@common.kspec_params
@click.option(
	'-o', '--output',
	required=True,
	type=common.filepath(writable=True),
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
@click.option('--progress/--no-progress', default=True, help="Show/don't show progress meter.")
@common.progress_arg()
@click.option('--dump-params', is_flag=True, hidden=True)
@click.pass_context
def create(ctx: click.Context,
           list_file: Optional[TextIO],
           ldir: Optional[str],
           files_arg: List[str],
           output: str,
           prefix: Optional[str],
           k: Optional[int],
           meta_file: Optional[TextIO],
           ids_file: Optional[TextIO],
           db_params: bool,
           progress: bool,
           dump_params: bool,
           ):
	"""Create k-mer signatures from genome sequences."""

	common.check_params_group(ctx, ['list_file', 'files_arg'], True, True)

	# Get sequence files
	ids, files = common.get_sequence_files(files_arg, list_file, ldir)
	seqfiles = SequenceFile.from_paths(files, 'fasta', 'auto')

	# Get kmerspec
	kspec = common.kspec_from_params(k, prefix)

	if db_params:
		if kspec is None:
			ctxobj = ctx.obj  # type: common.CLIContext
			ctxobj.require_signatures()
			kspec = ctx.obj.signatures.kmerspec
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
	sigs = calc_file_signatures(kspec, seqfiles, progress='click' if progress else None)
	sigs = AnnotatedSignatures(sigs, ids, meta)
	dump_signatures(output, sigs)
