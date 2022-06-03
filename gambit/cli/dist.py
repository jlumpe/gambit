import sys
from typing import Optional, TextIO, List

import click

from . import common
from .root import cli
from gambit.seq import SequenceFile
from gambit.sigs import load_signatures
from gambit.sigs.calc import calc_file_signatures
from gambit.metric import jaccarddist_matrix
import gambit.util.json as gjson
from gambit.util.progress import progress_config
from gambit.cluster import dump_dmat_csv


def fmt_kspec(kspec):
	return f'{kspec.k}/{kspec.prefix_str}'


@cli.command(name='dist', no_args_is_help=True)
@common.kspec_params
@click.option('-o', 'output', type=common.filepath(writable=True), required=True, help='Output file.')
@click.option('-q', type=common.filepath(exists=True), multiple=True, help='Query genome(s) (may be used multiple times).')
@click.option('--ql', type=click.File('r'), help='File containing paths to query genomes, one per line.')
@click.option('--qdir', type=common.dirpath(exists=True), default='.', help='Parent directory of files in --ql.')
@click.option('--qs', type=common.filepath(exists=True), help='Query signature file.')
@click.option('-r', type=common.filepath(exists=True), multiple=True, help='Reference genome (may be used multiple times).')
@click.option('--rl', type=click.File('r'), help='File containing paths to reference genomes, one per line.')
@click.option('--rdir', type=common.dirpath(exists=True), default='.', help='Parent directory of files in --rl.')
@click.option('--rs', type=common.filepath(exists=True), help='Reference signature file.')
@click.option('-s', '--square', is_flag=True, help='Calculate square distance matrix using query signatures only.')
@click.option('-d', '--use-db', is_flag=True, help='Use reference signatures from database.')
@common.progress_arg()
@click.option('--dump-params', is_flag=True, hidden=True)
@click.pass_context
def dist_cmd(ctx: click.Context,
             k: Optional[int],
             prefix: Optional[str],
             output: str,
             q: List[str],
             ql: Optional[TextIO],
             qdir: Optional[str],
             qs: Optional[str],
             r: List[str],
             rl: Optional[TextIO],
             rdir: Optional[str],
             rs: Optional[str],
             square: bool,
             use_db: bool,
             progress: bool,
             dump_params: bool,
             ):
	"""Calculate the GAMBIT distances between a set of query geneomes and a set of reference genomes.

	"""
	common.check_params_group(ctx, ['q', 'ql', 'qs'], True, True)
	common.check_params_group(ctx, ['r', 'rl', 'rs', 'use_db', 'square'], True, True)

	# Query files/signatures
	if qs is not None:
		query_sigs = load_signatures(qs)
		query_ids = query_sigs.ids
		query_files = None
	else:
		query_ids, query_files = common.get_sequence_files(q, ql, qdir)
		query_sigs = None

	common.warn_duplicate_file_ids(query_ids, 'Warning: the following query file IDs are present more than once: {ids}')

	# Ref files/signatures
	if rs is not None:
		ref_sigs = load_signatures(rs)
		ref_ids = ref_sigs.ids
		ref_files = None
	elif use_db:
		ctxobj = ctx.obj  # type: common.CLIContext
		ctxobj.require_signatures()
		ref_sigs = ctxobj.signatures
		ref_ids = ref_sigs.ids
		ref_files = None
	elif square:
		ref_ids = query_ids
		ref_files = ref_sigs = None
	else:
		ref_ids, ref_files = common.get_sequence_files(r, rl, rdir)
		ref_sigs = None

	common.warn_duplicate_file_ids(ref_ids, 'Warning: the following reference file IDs are present more than once: {ids}')

	# Kmerspec
	kspec = common.kspec_from_params(k, prefix)

	if kspec is None:
		if query_sigs is not None and ref_sigs is not None and query_sigs.kmerspec != ref_sigs.kmerspec:
			raise click.ClickException(
				f'K-mer search parameters of query signatures ({fmt_kspec(query_sigs.kmerspec)}) do '
				f'not match those of reference signatures ({fmt_kspec(ref_sigs.kmerspec)}).'
			)
		if query_sigs is not None:
			kspec = query_sigs.kmerspec
		elif ref_sigs is not None:
			kspec = ref_sigs.kmerspec
		else:
			raise click.ClickException('Require values for the -k and --prefix options.')

	else:
		if query_sigs is not None and query_sigs.kmerspec != kspec:
			raise click.ClickException(
				f'K-mer search parameters from command line options ({fmt_kspec(kspec)}) do not'
				f'match those of query signatures ({fmt_kspec(query_sigs.kmerspec)}).')
		if ref_sigs is not None and ref_sigs.kmerspec != kspec:
			raise click.ClickException(
				f'K-mer search parameters from command line options ({fmt_kspec(kspec)}) do not'
				f'match those of reference signatures ({fmt_kspec(ref_sigs.kmerspec)}).')

	prog = 'click' if progress else None

	# Dump parsed parameters
	if dump_params:
		params = dict(
			query_files=[f.path for f in query_files],
			query_sigs_file=qs,
			query_ids=query_ids,
			ref_files=[f.path for f in ref_files],
			ref_sigs_file=rs,
			ref_ids=ref_ids,
			kmerspec=kspec,
			square=square,
		)
		gjson.dump(params, sys.stdout)
		return

	# Calculate signatures if needed
	if query_sigs is None:
		query_sigfiles = SequenceFile.from_paths(query_files, 'fasta', 'auto')
		query_pconf = progress_config(prog, desc='Calculating query genome signatures') if len(query_files) > 1 else None
		query_sigs = calc_file_signatures(kspec, query_sigfiles, progress=query_pconf)

	dist_pconf = progress_config(prog, desc='Calculating distances')

	if square:
		dmat = jaccarddist_matrix(query_sigs, query_sigs, progress=dist_pconf)  # TODO use jaccarddist_pairwise

	else:
		if ref_sigs is None:
			ref_sigfiles = SequenceFile.from_paths(ref_files, 'fasta', 'auto')
			ref_pconf = progress_config('click', desc='Calculating reference genome signatures') if len(ref_files) > 1 else None
			ref_sigs = calc_file_signatures(kspec, ref_sigfiles, progress=ref_pconf)

		dmat = jaccarddist_matrix(query_sigs, ref_sigs, progress=dist_pconf)

	# Output
	dump_dmat_csv(output, dmat, query_ids, ref_ids)  # TODO different output formats
