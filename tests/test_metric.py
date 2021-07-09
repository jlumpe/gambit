"""Test distance metric calculations."""

import pickle

import pytest
import numpy as np

from gambit.metric import jaccard_sparse, jaccarddist_sparse, jaccard_bits, \
	jaccard_generic, jaccard_sparse_array, SCORE_DTYPE, BOUNDS_DTYPE
from gambit.kmers import sparse_to_dense
from gambit.signatures import SignatureArray
from gambit.test import make_signatures


@pytest.fixture(scope='module')
def load_test_coords_col(test_data):
	"""Function to load k-mer coordinates from file."""

	def load_test_coords_col_func():
		with open(test_data / 'kmer_coords/coords.pickle', 'rb') as fobj:
			signatures_list = pickle.load(fobj)

		return SignatureArray(signatures_list)

	return load_test_coords_col_func


@pytest.fixture(params=[
	(4, 'u2'),
	(4, 'i2'),
	(7, 'u2'),
	(7, 'i2'),
	(9, 'u4'),
	(9, 'i4'),
	(9, 'u8'),
	(9, 'i8'),
	None,
])
def coords_params(request, load_test_coords_col):
	"""Tuple of (k, SignatureArray) to test on."""

	if request.param is None:
		# Load coords from file
		k = 11
		sigs = load_test_coords_col()

	else:
		# Create coords
		k, dtype = request.param

		np.random.seed(0)
		sigs = make_signatures(k, 25, dtype)

	return k, sigs


def test_jaccard_single(coords_params):
	"""Test calculating single scores at a time."""

	k, sigs = coords_params

	# Iterate over all pairs
	for i, coords1 in enumerate(sigs):
		vec1 = sparse_to_dense(k, coords1)
		for j, coords2 in enumerate(sigs):
			vec2 = sparse_to_dense(k, coords2)

			score = jaccard_sparse(coords1, coords2)

			# Check range
			assert 0 <= score <= 1

			# Check vs slow version
			assert np.isclose(score, jaccard_generic(coords1, coords2))

			# Check distance
			assert np.isclose(jaccarddist_sparse(coords1, coords2), 1 - score)

			# Check dense bit vector version
			assert np.isclose(jaccard_bits(vec1, vec2), score)

			# Check score vs. self is one (unless empty)
			if i == j:
				assert score == 0 if len(coords1) == 0 else 1


@pytest.mark.parametrize('alt_bounds_dtype', [False, True])
def test_jaccard_sparse_array(coords_params, alt_bounds_dtype):
	"""Test jaccard_sparse_array() function."""

	k, sigs = coords_params

	# The inner Cython function takes a specific type for the bounds array.
	# Try with this type and a different type, should be converted automatically by the outer Python func
	if alt_bounds_dtype:
		sigs = SignatureArray.from_arrays(sigs.values, sigs.bounds.astype('i4'))
		assert sigs.bounds.dtype != BOUNDS_DTYPE
	else:
		assert sigs.bounds.dtype == BOUNDS_DTYPE

	for i, coords1 in enumerate(sigs):
		scores = jaccard_sparse_array(coords1, sigs)
		assert scores.shape == (len(sigs),)

		# Check against single coords
		for j, coords2 in enumerate(sigs):
			assert scores[j] == jaccard_sparse(coords1, coords2)

		# Check distance
		dists = jaccard_sparse_array(coords1, sigs, distance=True)
		assert np.allclose(dists, 1 - scores)

	# Check pre-allocated output
	out = np.empty(len(sigs), dtype=SCORE_DTYPE)
	jaccard_sparse_array(sigs[0], sigs, out=out)
	assert np.array_equal(out, jaccard_sparse_array(sigs[0], sigs))

	# Wrong size
	out2 = np.empty(len(sigs) + 1, dtype=SCORE_DTYPE)
	with pytest.raises(ValueError):
		jaccard_sparse_array(sigs[0], sigs, out=out2)

	# Wrong dtype
	out3 = np.empty(len(sigs), dtype=int)
	with pytest.raises(ValueError):
		jaccard_sparse_array(sigs[0], sigs, out3)


def test_different_dtypes():
	"""Test metric on sparse arrays with different dtypes."""

	np.random.seed(0)
	dtypes = ['i2', 'u2', 'i4', 'u4', 'i8', 'u8']

	# Test all pairs of dtypes
	for i, dt1 in enumerate(dtypes):
		for j in range(i + 1, len(dtypes)):
			dt2 = dtypes[j]

			# Bit length of largest positive integer both dtypes can store
			nbits = min(np.iinfo(dt).max.bit_length() for dt in [dt1, dt2])

			# Try for k == 8, but not larger than will fit in both dtypes
			k = min(nbits // 4, 8)

			# Create test signatures as u8 dtype (will fit everything)
			sigs = make_signatures(k, 5, 'u8')

			# Convert signatures to each dtype, making sure there is no overflow
			try:
				# Tell numpy to raise error on overflow
				old_err = np.seterr(over='raise')

				sigs1 = SignatureArray(sigs, dtype=dt1)
				sigs2 = SignatureArray(sigs, dtype=dt2)

			finally:
				np.seterr(**old_err)

			assert sigs1.values.dtype == dt1
			assert sigs2.values.dtype == dt2

			for k in range(len(sigs)):
				# Test individually
				for l in range(k + 1, len(sigs)):
					expected = jaccard_sparse(sigs[k], sigs[l])
					assert jaccard_sparse(sigs1[k], sigs2[l]) == expected
					assert jaccard_sparse(sigs2[k], sigs1[l]) == expected

				# Test against all of the other using jaccard_sparse_array
				expected_all = jaccard_sparse_array(sigs[k], sigs)
				assert np.array_equal(jaccard_sparse_array(sigs1[k], sigs2), expected_all)
				assert np.array_equal(jaccard_sparse_array(sigs2[k], sigs1), expected_all)
