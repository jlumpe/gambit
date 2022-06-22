"""Distance matrices and basic clustering/trees."""

from typing import Union, Optional, Sequence, TextIO, Tuple, List
import csv

import numpy as np
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage
from Bio.Phylo.BaseTree import Tree, Clade

from gambit.util.io import FilePath, maybe_open
from gambit.util.misc import zip_strict


def hclust(dmat: np.ndarray) -> np.ndarray:
	"""Perform hierarchical clustering on a distance matrix.

	Uses the UPGMA linkage method.

	Parameters
	----------
	dmat
		Pairwise distance matrix.

	Returns
	-------
	np.ndarray
		Linkage matrix as returned by :func:`scipy.cluster.hierarchy.linkage`.
	"""
	assert dmat.ndim == 2
	sm = squareform(dmat)
	return linkage(sm, method='average')


def linkage_to_bio_tree(link: np.ndarray, labels: Sequence[str]) -> Tree:
	"""Convert SciPy linkage matrix to BioPython phylogenetic tree object.

	Parameters
	----------
	nodes
		Matrix returned by :func:`scipy.cluster.hierarchy.linkage`.
	labels
		Leaf node names.

	Returns
	-------
	Bio.Phylo.BaseTree.Tree
		BioPython tree object.
	"""
	nleaves = link.shape[0] + 1
	nnodes = nleaves * 2 + 1

	labels = list(labels)
	assert len(labels) == nleaves

	# Make leaves
	clades = [Clade(name=name) for name in labels]

	# Make internal nodes
	for left_i, right_i, height, size in link:
		left_i = int(left_i)
		left = clades[left_i]
		left.branch_length = height - (0 if left_i < nleaves else link[left_i - nleaves, 2])

		right_i = int(right_i)
		right = clades[right_i]
		right.branch_length = height - (0 if right_i < nleaves else link[right_i - nleaves, 2])

		clades.append(Clade(clades=[left, right]))

	return Tree(root=clades[-1], rooted=True)


def check_tree_matches_linkage(tree: Tree, link: np.ndarray, labels, atol=1e-5):
	"""
	Testing function to check that the output of :func:`.linkage_to_bio_tree` is consistent with its
	input.

	Raises
	------
	AssertionError
		If any checks fail.
	"""
	nleaves = len(labels)
	nnodes = nleaves * 2 - 1
	# assert tree.rooted
	assert link.shape[0] == nleaves - 1

	label_to_index = {l: i for i, l in enumerate(labels)}
	child_ids = [{int(l), int(r)} for l, r, h, s in link]

	def height_close(h1, h2): return np.isclose(h1, h2, rtol=0, atol=atol)

	def check_clade(clade):
		if clade.clades:  # Internal node
			left, right = clade.clades
			left_i, left_h = check_clade(left)
			right_i, right_h = check_clade(right)

			# Find node ID/index based on IDs of children
			try:
				i = child_ids.index({left_i, right_i}) + nleaves
			except ValueError:
				raise AssertionError(f'Linkage matrix contains no node with child IDs {left_i} and {right_i}') from None

			assert int(link[i - nleaves, 0]) == left_i
			assert int(link[i - nleaves, 1]) == right_i
			height = link[i - nleaves, 2]
			assert height_close(left_h + left.branch_length, height)
			assert height_close(right_h + right.branch_length, height)

		else:  # Leaf node
			i = label_to_index[clade.name]
			height = 0

		return i, height

	root_i, root_h = check_clade(tree.root)
	assert root_i == nleaves * 2 - 2


def dump_dmat_csv(file: Union[FilePath, TextIO],
                  dmat: np.ndarray,
                  row_ids: Sequence,
                  col_ids: Sequence,
                  corner: Optional[str] = None,
                  fmt: str = '0.4f',
                  ):
	"""Write distance matrix to file in CSV format."""

	with maybe_open(file, 'w', newline='') as fobj:
		writer = csv.writer(fobj)
		writer.writerow([corner or '', *map(str, col_ids)])
		for row_id, values in zip_strict(row_ids, dmat):
			values_str = (format(d, fmt) for d in values)
			writer.writerow([str(row_id), *values_str])


def load_dmat_csv(file: Union[FilePath, TextIO]) -> Tuple[np.ndarray, List[str], List[str]]:
	"""Load distance matrix from CSV file.

	Returns
	-------
	Tuple
		``(matrix, row_ids, col_ids)`` tuple.
	"""

	with maybe_open(file, newline='') as fobj:
		reader = csv.reader(fobj)
		col_ids = next(reader)[1:]
		nc = len(col_ids)

		rows = []
		row_ids = []

		for rid, *row_str in reader:
			row_ids.append(rid)
			rows.append(np.fromiter(map(float, row_str), np.float32, count=nc))

		values = np.stack(rows)
		return values, row_ids, col_ids
