import sys
from typing import Optional, TextIO, List

import click
from Bio import Phylo

from . import common
from .root import cli
from gambit.seq import SequenceFile
from gambit.sigs import load_signatures
from gambit.sigs.calc import calc_file_signatures
from gambit.metric import jaccarddist_pairwise
from gambit.util.progress import progress_config
from gambit._cython.threads import omp_set_num_threads
from gambit.cluster import hclust, linkage_to_bio_tree


@cli.command(name='tree', no_args_is_help=True)
@common.genome_files_arg()
@click.option(
	'-l', 'listfile',
	type=click.File('r'),
	metavar='LISTFILE',
	help='File containing paths to genomes.',
)
@click.option('--ldir', type=common.dirpath(), default='.', help='Parent directory of paths in LISTFILE.')
@click.option(
	'-s', '--sigfile',
	type=common.filepath(exists=True),
	help='GAMBIT signatures file.',
)
@common.kspec_params
@click.option('-c', '--cores', type=click.IntRange(min=1), help='Number of CPU cores to use.')
@common.progress_arg()
@click.pass_context
def tree_cmd(ctx: click.Context,
             listfile: Optional[TextIO],
             ldir: Optional[str],
             files_arg: List[str],
             sigfile: Optional[str],
             k: Optional[int],
             prefix: Optional[str],
             progress: bool,
             cores: Optional[int],
             ):
	"""
	Estimate a relatedness tree for a set of genomes and output in Newick format.
	"""
	common.check_params_group(ctx, ['files_arg', 'listfile', 'sigfile'], True, True)

	if cores is not None:
		omp_set_num_threads(cores)

	pconf = progress_config('click' if progress else None)

	# Files/signatures
	if sigfile is not None:
		sigs = load_signatures(sigfile)
		labels = sigs.ids

	else:
		labels, genome_files = common.get_sequence_files(files_arg, listfile, ldir)
		common.warn_duplicate_file_ids(labels, 'Warning: the following query file IDs are present more than once: {ids}')

		kspec = common.kspec_from_params(k, prefix)
		if kspec is None:
			raise click.ClickException('Require values for the -k and --prefix options.')

		sigfiles = SequenceFile.from_paths(genome_files, 'fasta', 'auto')
		sigs = calc_file_signatures(kspec, sigfiles, progress=pconf.update(desc='Calculating signatures'), max_workers=cores)

	# Calculate distances
	dmat = jaccarddist_pairwise(sigs, progress=pconf.update(desc='Calculating distances'))

	# Cluster
	link = hclust(dmat)
	tree = linkage_to_bio_tree(link, labels)

	Phylo.write(tree, sys.stdout, 'newick')
