"""Test the gambit.cluster module."""

from gambit import cluster

import numpy as np


def test_dmat_csv(tmp_path):
	"""Test dumping and loading distance matrix in CSV format."""

	nr = 10
	nc = 20
	row_ids = [f'r{i}' for i in range(nr)]
	col_ids = [f'c{i}' for i in range(nc)]
	dmat = np.random.rand(nr, nc)

	file = tmp_path / 'dmat.csv'
	cluster.dump_dmat_csv(file, dmat, row_ids, col_ids)
	dmat2, rids2, cids2 = cluster.load_dmat_csv(file)

	assert np.allclose(dmat, dmat2, atol=1e-4)
	assert rids2 == row_ids
	assert cids2 == col_ids


def test_tree():
	"""Test hierarchical clustering and converting to BioPython tree object."""

	# Made-up linkage array, (((A, B), C), D)
	link = np.asarray([
		[0, 1, .1, 2],
		[2, 4, .2, 3],
		[3, 5, .3, 4],
	], dtype=float)
	labels = list('ABCD')

	tree = cluster.linkage_to_bio_tree(link, labels)
	assert tree.rooted
	cluster.check_tree_matches_linkage(tree, link, labels)
