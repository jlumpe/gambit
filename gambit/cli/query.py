import sys
from typing import List, TextIO

import click

from .context import CLIContext
from gambit.db import GAMBITDatabase
from gambit.query import QueryParams, query_parse
from gambit.io.seq import SequenceFile
from gambit.util.progress import ClickProgressMeter


def get_exporter(outfmt: str):
	if outfmt == 'csv':
		from gambit.io.export.csv import CSVResultsExporter
		return CSVResultsExporter()

	if outfmt == 'json':
		from gambit.io.export.json import JSONResultsExporter
		return JSONResultsExporter()

	if outfmt == 'archive':
		from gambit.io.export.archive import ResultsArchiveWriter
		return ResultsArchiveWriter(install_info=True)

	assert 0


@click.command()
@click.option(
	'-o', '--output',
	type=click.File(mode='w'),
	default=sys.stdout,
	help='File path to write to. If omitted will write to stdout.',
)
@click.option(
	'-s', '--seqfmt',
	type=click.Choice(['fasta']),
	default='fasta',
	help='Format of sequence files. Currently only FASTA is supported.',
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
@click.argument(
	'files',
	nargs=-1,
	type=click.Path(exists=True, dir_okay=False),
	required=True,
	metavar='GENOMES...',
)
@click.pass_obj
def query(ctxobj: CLIContext, files: List[str], output: TextIO, seqfmt: str, outfmt: str, strict: bool):
	"""Predict taxonomy of microbial samples from genome sequences."""
	gset = ctxobj.genomeset()
	ref_sigs = ctxobj.signatures()
	db = GAMBITDatabase(gset, ref_sigs)

	params = QueryParams(classify_strict=strict)
	files = SequenceFile.from_paths(files, seqfmt)
	exporter = get_exporter(outfmt)

	results = query_parse(db, files, params, progress=ClickProgressMeter)
	exporter.export(output, results)
