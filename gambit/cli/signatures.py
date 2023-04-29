from typing import Optional, TextIO, List
import sys

import click

from . import common
from .root import cli
import gambit.util.json as gjson
from gambit.sigs import SignaturesMeta, AnnotatedSignatures, load_signatures, dump_signatures
from gambit.sigs.calc import calc_file_signatures
from gambit.util.io import read_lines
from gambit.kmers import DEFAULT_KMERSPEC


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
			hdf5_format_version=sigs.format_version,  # HDF5-specific
		)
		kw = dict(indent=' ', sort_keys=True) if pretty else dict()
		gjson.dump(data, sys.stdout, **kw)

	else:
		rows1 = [
			('Genome Count:', len(sigs)),
			('k:', sigs.kmerspec.k),
			('Prefix:', sigs.kmerspec.prefix.decode('ascii')),
			('File format:', f'HDF5, version {sigs.format_version}'),  # HDF5-specific
			('Data type:', sigs.dtype),
		]
		common.print_table(rows1, colsep='  ')

		print('Metadata:')

		rows2 = [
			('ID:', format_none(sigs.meta.id)),
			('Name:', format_none(sigs.meta.name)),
			('Version:', format_none(sigs.meta.version)),
			('Description:', format_none(sigs.meta.description)),
			('Genome ID attribute:', format_none(sigs.meta.id_attr)),
			('Has extra:', 'yes' if sigs.meta.extra else 'no'),
		]
		common.print_table(rows2, colsep='  ', left='  ')


@signatures_group.command(no_args_is_help=True)
@common.genome_files_arg()
@common.listfile_param('-l', 'listfile', metavar='LISTFILE', help='File containing paths to genome files, one per line.')
@common.listfile_dir_param('--ldir', file_metavar='LISTFILE')
@common.kspec_params()
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
@common.progress_param()
@common.cores_param()
@click.option('--dump-params', is_flag=True, hidden=True)
@click.pass_context
def create(ctx: click.Context,
           listfile: Optional[TextIO],
           ldir: Optional[str],
           files_arg: List[str],
           output: str,
           prefix: Optional[str],
           k: Optional[int],
           meta_file: Optional[TextIO],
           ids_file: Optional[TextIO],
           db_params: bool,
           progress: bool,
           cores: Optional[int],
           dump_params: bool,
           ):
	"""Create k-mer signatures from genome sequences."""

	common.check_params_group(ctx, ['listfile', 'files_arg'], True, True)

	# Get sequence files
	ids, files = common.get_sequence_files(files_arg, listfile, ldir)

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
		kspec = DEFAULT_KMERSPEC

	# Metadata / IDs
	if meta_file is not None:
		meta = gjson.load(meta_file, SignaturesMeta)
	else:
		meta = None

	if ids_file is not None:
		ids = list(read_lines(ids_file))
		if len(ids) != len(files):
			raise click.ClickException(f'Number of IDs ({len(ids)}) does not match number of genomes ({len(files)}).')

	common.warn_duplicate_file_ids(ids, 'Warning: the following file IDs are present more than once: {ids}')

	# Dump parameters
	if dump_params:
		params = dict(
			kmerspec=kspec,
			files=[f.path for f in files],
			meta=meta,
			ids=ids,
		)
		gjson.dump(params, sys.stdout)
		return

	# Calculate and save
	sigs = calc_file_signatures(kspec, files, progress='click' if progress else None, max_workers=cores)
	sigs = AnnotatedSignatures(sigs, ids, meta)
	dump_signatures(output, sigs)
