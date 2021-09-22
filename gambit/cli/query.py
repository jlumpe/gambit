import sys
from typing import TextIO

import click

from .common import CLIContext, seq_file_params, get_seq_files
from gambit.db import GAMBITDatabase
from gambit.query import QueryParams, query_parse
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


@click.command(name='query')
@seq_file_params()
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
@click.pass_obj
def query_cmd(ctxobj: CLIContext,
              output: TextIO,
              outfmt: str,
              strict: bool,
              **kw,
              ):
	"""Predict taxonomy of microbial samples from genome sequences."""
	gset = ctxobj.genomeset()
	ref_sigs = ctxobj.signatures()
	db = GAMBITDatabase(gset, ref_sigs)

	seqfiles = get_seq_files(kw)
	params = QueryParams(classify_strict=strict)
	exporter = get_exporter(outfmt)

	results = query_parse(db, seqfiles, params, progress=ClickProgressMeter)
	exporter.export(output, results)
