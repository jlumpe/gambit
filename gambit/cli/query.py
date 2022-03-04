import sys
from typing import TextIO, Optional, List

import click

from .common import CLIContext, genome_files_arg, filepath
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


@click.command(name='query')
@genome_files_arg()
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
	'--sigfile',
	type=filepath(exists=True),
	help='File containing query signatures, to use in place of GENOMES.',
)
@click.pass_obj
def query_cmd(ctxobj: CLIContext,
              files: List[str],
              sigfile: Optional[str],
              output: TextIO,
              outfmt: str,
              strict: bool,
              ):
	"""Predict taxonomy of microbial samples from genome sequences."""

	if sigfile and files:
		raise click.ClickException('The --sigfile option is mutually exclusive with GENOMES')
	if not sigfile and not files:
		raise click.ClickException('Must supply at least one genome file or a value for --sigfile.')

	db = ctxobj.get_database()
	seqfiles = SequenceFile.from_paths(files, 'fasta', 'auto')
	params = QueryParams(classify_strict=strict)
	exporter = get_exporter(outfmt)
	pconf = progress_config('click', file=sys.stderr)

	if sigfile:
		sigs = load_signatures(sigfile)
		inputs = [QueryInput(id) for id in sigs.ids]
		results = query(db, sigs, params, inputs=inputs, progress=pconf)

	else:
		results = query_parse(db, seqfiles, params, progress=pconf)

	exporter.export(output, results)
