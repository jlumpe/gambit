"""Test distance metric calculations."""

import pickle

import pytest
import numpy as np

from gambit.metric import jaccard_sparse, jaccarddist_sparse, jaccard_bits, \
	jaccard_generic, jaccard_sparse_array, jaccard_sparse_matrix, SCORE_DTYPE, BOUNDS_DTYPE
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
	(7, 'u2'),
	(9, 'u4'),
	(9, 'i4'),
	(9, 'u8'),
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
		sigs = make_signatures(k, 40, dtype)

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


class TestJaccardSparseMatrix:
	"""Test the jaccard_sparse_matrix() function."""

	def make_output_array(self, refs, queries, dtype=SCORE_DTYPE):
		return np.empty((len(refs), len(queries)), dtype=dtype)

	@pytest.fixture()
	def queries(self, coords_params):
		"""queries argument."""
		k, sigs = coords_params
		# Make it a little different than the reference signatures
		return sigs[::2]

	@pytest.fixture()
	def refs_array(self, coords_params):
		"""refs argument as SignatureArray."""
		k, sigs = coords_params
		# Make it a little different than the query signatures
		return sigs[5:]

	@pytest.fixture()
	def refs(self, refs_array):
		"""refs argument.
		TODO: test other AbstractSignaturesArray types
		"""
		return refs_array

	@pytest.fixture()
	def expected(self, queries, refs):
		"""Expected array of scores, calculated one at a time."""
		scores = np.empty((len(queries), len(refs)), dtype=SCORE_DTYPE)

		for i, q in enumerate(queries):
			for j, r in enumerate(refs):
				scores[i, j] = jaccard_sparse(queries[i], refs[j])

		return scores

	@pytest.mark.parametrize('use_ref_indices', [False, True])
	@pytest.mark.parametrize('chunksize', [None, 10])
	def test_basic(self, queries, refs, expected, use_ref_indices, chunksize):

		if use_ref_indices:
			ref_indices = [i for i in range(len(refs)) if i % 3 != 0]
			expected = expected[:, ref_indices]
		else:
			ref_indices = None

		scores = jaccard_sparse_matrix(queries, refs, ref_indices=ref_indices, chunksize=chunksize)
		assert np.array_equal(scores, expected)

	def test_distance(self, queries, refs, expected):
		dists = jaccard_sparse_matrix(queries, refs, distance=True)
		assert np.allclose(dists, 1 - expected)

	def test_out(self, queries, refs, expected):
		"""Test using pre-allocated output array."""
		out = np.empty((len(queries), len(refs)), dtype=SCORE_DTYPE)
		jaccard_sparse_matrix(queries, refs, out=out)
		assert np.array_equal(out, expected)

		# Wrong size
		out2 = np.empty((len(refs) + 1, len(queries)), dtype=SCORE_DTYPE)
		with pytest.raises(ValueError):
			jaccard_sparse_matrix(queries, refs, out=out2)

		# Wrong dtype
		out3 = np.empty((len(refs) + 1, len(queries)), dtype=int)
		with pytest.raises(ValueError):
			jaccard_sparse_array(queries, refs, out3)

	def test_progress(self, queries, refs, expected):
		pass  # TODO


def test_different_dtypes():
	"""Test metric on sparse arrays with different dtypes."""

	np.random.seed(0)
	dtypes = ['u2', 'i4', 'u4', 'u8']

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

			# Expected all-by-all distance matrix from u8 signatures
			expected = jaccard_sparse_matrix(sigs, sigs)

			for k in range(len(sigs)):
				# Test individually
				for l in range(k + 1, len(sigs)):
					assert jaccard_sparse(sigs1[k], sigs2[l]) == expected[k, l]
					assert jaccard_sparse(sigs2[k], sigs1[l]) == expected[k, l]

				# Test each against all of the others using jaccard_sparse_array
				assert np.array_equal(jaccard_sparse_array(sigs1[k], sigs2), expected[k, :])
				assert np.array_equal(jaccard_sparse_array(sigs2[k], sigs1), expected[k, :])

			# Test full matrix using jaccard_sparse_matrix
			assert np.array_equal(jaccard_sparse_matrix(sigs1, sigs2), expected)
