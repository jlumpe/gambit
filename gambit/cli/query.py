import sys
from typing import TextIO, Optional, List

import click

from . import common
from .root import cli
from gambit.query import QueryParams, QueryInput, query, query_parse
from gambit.util.progress import progress_config
from gambit.sigs import load_signatures
from gambit._cython.threads import omp_set_num_threads


def get_exporter(outfmt: str):
	if outfmt == 'csv':
		from gambit.results.csv import CSVResultsExporter
		return CSVResultsExporter()

	if outfmt == 'json':
		from gambit.results.json import JSONResultsExporter
		return JSONResultsExporter()

	if outfmt == 'archive':
		from gambit.results.archive import ResultsArchiveWriter
		return ResultsArchiveWriter(install_info=True)

	assert 0


@cli.command(name='query', no_args_is_help=True)
@common.genome_files_arg()
@click.option(
	'-l', 'listfile',
	type=click.File('r'),
	metavar='LISTFILE',
	help='File containing paths to genomes.',
)
@click.option('--ldir', type=common.dirpath(), default='.', help='Parent directory of paths in LISTFILE.')
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
@common.progress_arg()
@click.option('-c', '--cores', type=click.IntRange(min=1), help='Number of CPU cores to use.')
@click.pass_context
def query_cmd(ctx: click.Context,
              listfile: Optional[TextIO],
              ldir: Optional[str],
              files_arg: List[str],
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
		inputs = [QueryInput(id) for id in sigs.ids]
		results = query(db, sigs, params, inputs=inputs, progress=pconf)

	else:
		ids, files = common.get_sequence_files(files_arg, listfile, ldir)
		common.warn_duplicate_file_ids(ids, 'Warning: the following query file IDs are present more than once: {ids}')
		results = query_parse(
			db, files, params,
			file_labels=ids,
			progress=pconf,
			parse_kw=dict(max_workers=cores),
		)

	exporter.export(output, results)
