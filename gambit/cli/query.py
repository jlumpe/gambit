import sys
from typing import TextIO, Optional

import click

from .common import CLIContext, seq_file_params, get_seq_files
from gambit.query import QueryParams, QueryInput, query, query_parse
from gambit.util.progress import ClickProgressMeter
from gambit.signatures.hdf5 import HDF5Signatures


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
@click.option(
	'--sigfile',
	type=click.Path(exists=True, dir_okay=False, readable=True),
	help='File containing query signatures, to use in place of GENOMES.',
)
@click.pass_obj
def query_cmd(ctxobj: CLIContext,
              sigfile: Optional[str],
              output: TextIO,
              outfmt: str,
              strict: bool,
              **kw,
              ):
	"""Predict taxonomy of microbial samples from genome sequences."""

	db = ctxobj.get_database()
	seqfiles = get_seq_files(kw)
	params = QueryParams(classify_strict=strict)
	exporter = get_exporter(outfmt)

	if sigfile and seqfiles:
		raise click.ClickException('The --sigfile option is mutually exclusive with GENOMES')

	elif sigfile:
		with HDF5Signatures.open(sigfile) as sigfile:
			sigs = sigfile[:]
		inputs = [QueryInput(id) for id in sigfile.ids]
		results = query(db, sigs, params, inputs=inputs, progress=ClickProgressMeter)

	elif seqfiles:
		results = query_parse(db, seqfiles, params, progress=ClickProgressMeter)

	else:
		raise click.ClickException('Must supply at least one genome file or a value for --sigfile.')

	exporter.export(output, results)
