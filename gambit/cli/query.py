import sys
from typing import TextIO, Optional, List

import click

from .common import CLIContext, genome_files_arg, filepath, dirpath, read_genomes_list_file
from .root import cli
from gambit.query import QueryParams, QueryInput, query, query_parse
from gambit.util.progress import progress_config
from gambit.sigs import load_signatures
from gambit.seq import SequenceFile


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
@genome_files_arg()
@click.option(
	'-l', 'listfile',
	type=click.File('r'),
	metavar='LISTFILE',
	help='File containing paths to genomes.',
)
@click.option('--ldir', type=dirpath(), default='.', help='Parent directory of paths in LISTFILE.')
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
	type=filepath(exists=True),
	help='File containing query signatures, to use in place of GENOMES.',
)
@click.pass_obj
def query_cmd(ctxobj: CLIContext,
              listfile: Optional[TextIO],
              ldir: Optional[str],
              files: List[str],
              sigfile: Optional[str],
              output: TextIO,
              outfmt: str,
              strict: bool,
              ):
	"""Predict taxonomy of microbial samples from genome sequences."""

	count = (len(files) > 0) + (listfile is not None) + (sigfile is not None)
	if count == 0:
		raise click.ClickException('Must give value(s) for one of GENOMES, -s, or -l.')
	if count > 1:
		raise click.ClickException('The GENOMES, -s, and -l parameters are mutally exclusive.')

	db = ctxobj.get_database()
	params = QueryParams(classify_strict=strict)
	exporter = get_exporter(outfmt)
	pconf = progress_config('click', file=sys.stderr)

	if sigfile:
		sigs = load_signatures(sigfile)
		inputs = [QueryInput(id) for id in sigs.ids]
		results = query(db, sigs, params, inputs=inputs, progress=pconf)

	else:
		if listfile is not None:
			files = read_genomes_list_file(listfile, ldir)
		seqfiles = SequenceFile.from_paths(files, 'fasta', 'auto')
		results = query_parse(db, seqfiles, params, progress=pconf)

	exporter.export(output, results)
