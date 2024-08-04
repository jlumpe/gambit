import sys
from typing import TextIO, Optional

import click

from . import common
from .root import cli
from gambit.query import QueryParams, query, query_parse
from gambit.util.progress import progress_config
from gambit.sigs import load_signatures
from gambit.results import CSVResultsExporter, JSONResultsExporter, ResultsArchiveWriter
from gambit._cython.threads import omp_set_num_threads


def get_exporter(outfmt: str):
	if outfmt == 'csv':
		return CSVResultsExporter()

	if outfmt == 'json':
		return JSONResultsExporter()

	if outfmt == 'archive':
		return ResultsArchiveWriter()

	raise ValueError(f'Invalid output format: {outfmt!r}')


@cli.command(name='query', no_args_is_help=True)
@common.genome_files_arg()
@common.listfile_param('-l', 'listfile', metavar='LISTFILE', help='File containing paths to query genomes, one per line.')
@common.listfile_dir_param('--ldir', file_metavar='LISTFILE')
@click.option(
	'-o', '--output',
	type=click.File(mode='w'),
	default=sys.stdout,
	help='File path to write to. If omitted will write to stdout.',
)
@click.option(
	'-f', '--outfmt',
	type=click.Choice(['csv', 'json', 'archive']),
	default='csv',
	help='Format to output results in.',
)
@click.option(
	'--strict/--no-strict',
	default=False,
	hidden=True,
)
@click.option(
	'-s', '--sigfile',
	type=common.filepath(exists=True),
	help='File containing query signatures, to use in place of GENOMES.',
)
@common.progress_param()
@common.cores_param()
@click.pass_context
def query_cmd(ctx: click.Context,
              listfile: Optional[TextIO],
              ldir: Optional[str],
              files_arg: list[str],
              sigfile: Optional[str],
              output: TextIO,
              outfmt: str,
              strict: bool,
              progress: bool,
              cores: Optional[int],
              ):
	"""Predict taxonomy of microbial samples from genome sequences."""

	common.check_params_group(ctx, ['files_arg', 'listfile', 'sigfile'], True, True)

	db = ctx.obj.get_database()
	params = QueryParams(classify_strict=strict)
	exporter = get_exporter(outfmt)
	pconf = progress_config('click' if progress else None)

	if cores is not None:
		omp_set_num_threads(cores)

	if sigfile:
		sigs = load_signatures(sigfile)
		results = query(db, sigs, params, progress=pconf)

	else:
		ids, files = common.get_sequence_files(files_arg, listfile, ldir)
		common.warn_duplicate_file_ids(ids, 'Warning: the following query file IDs are present more than once: {ids}')
		results = query_parse(
			db, files, params,
			labels=ids,
			progress=pconf,
			parse_kw=dict(max_workers=cores),
		)

	exporter.export(output, results)
