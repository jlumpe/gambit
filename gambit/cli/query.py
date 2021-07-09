import sys

import click

from .context import CLIContext
from gambit.db import GAMBITDatabase
from gambit.io.seq import SequenceFile, find_kmers_in_files
from gambit.query import runquery


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
@click.argument(
	'files',
	nargs=-1,
	type=click.Path(exists=True, dir_okay=False),
	required=True,
)
@click.pass_obj
def query(ctxobj: CLIContext, files, output, seqfmt: str, outfmt: str):
	"""Query database."""
	gset = ctxobj.genomeset()
	sigs = ctxobj.signatures()
	db = GAMBITDatabase(gset, sigs)

	# Parse files
	files = SequenceFile.from_paths(files, seqfmt)
	sigs = find_kmers_in_files(db.kmerspec, files)

	# Run query
	results = runquery(db, sigs, inputs=files)

	# Export results
	if outfmt == 'json':
		from gambit.query.export.json import JSONResultsExporter
		exporter = JSONResultsExporter()
	else:
		from gambit.query.export.csv import CSVResultsExporter
		exporter = CSVResultsExporter()

	exporter.export(output, results)
