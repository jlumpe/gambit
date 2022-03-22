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
