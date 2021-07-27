import sys

import click

from .context import CLIContext
from gambit.db import GAMBITDatabase
from gambit.query import QueryParams, query_parse
from gambit.io.seq import SequenceFile
from gambit.util.progress import ClickProgressMeter


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
	type=click.Choice(['csv', 'json']),
	default='json',
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
def query(ctxobj: CLIContext, files, output, seqfmt: str, outfmt: str, strict: bool):
	"""Predict taxonomy of microbial samples from genome sequences."""
	gset = ctxobj.genomeset()
	ref_sigs = ctxobj.signatures()
	db = GAMBITDatabase(gset, ref_sigs)

	params = QueryParams(classify_strict=strict)
	files = SequenceFile.from_paths(files, seqfmt)

	# Run query
	results = query_parse(db, files, params, progress=ClickProgressMeter)

	# Export results
	if outfmt == 'json':
		from gambit.export.json import JSONResultsExporter
		exporter = JSONResultsExporter()
	else:
		from gambit.export.csv import CSVResultsExporter
		exporter = CSVResultsExporter()

	exporter.export(output, results)
