"""Test the gambit.sigs.convert module."""

import pytest
import numpy as np

from gambit.sigs.convert import dense_to_sparse, sparse_to_dense, can_convert, \
	check_can_convert, convert_dense, convert_sparse
from gambit.kmers import KmerSpec
from gambit.test import random_seq


def test_dense_sparse_conversion():
	"""Test conversion between dense and sparse representations of k-mer coordinates."""

	for k in range(1, 10):

		kspec = KmerSpec(k, 'ATGAC')

		# Create dense signature with every 3rd k-mer
		vec = np.zeros(kspec.nkmers, dtype=bool)
		vec[np.arange(vec.size) % 3 == 0] = True

		# Convert to sparse
		sig = dense_to_sparse(vec)

		assert len(sig) == vec.sum()
		for index in sig:
			assert vec[index]

		# Check sorted
		assert np.all(np.diff(sig) > 0)

		# Check converting back
		assert np.array_equal(vec, sparse_to_dense(kspec, sig))


class TestKmerSpecConversion:
	"""Test converting signatures from one KmerSpec to another."""

	def test_can_convert(self):
		from_kspec = KmerSpec(11, 'ATGAC')

		compatible = [
			KmerSpec(11, 'ATGAC'),
			KmerSpec(8, 'ATGAC'),
			KmerSpec(10, 'ATGACA'),
			KmerSpec(8, 'ATGACA'),
		]

		for to_kspec in compatible:
			assert can_convert(from_kspec, to_kspec)
			check_can_convert(from_kspec, to_kspec)

		incompatible = [
			KmerSpec(11, 'CAGTA'),
			KmerSpec(12, 'ATGAC'),
			KmerSpec(11, 'ATGA'),
			KmerSpec(11, 'ATGACT'),
		]

		for to_kspec in incompatible:
			assert not can_convert(from_kspec, to_kspec)
			with pytest.raises(ValueError):
				check_can_convert(from_kspec, to_kspec)

	@pytest.fixture(scope='class')
	def seqs(self):
		np.random.seed(0)
		return [random_seq(100_000) for _ in range(100)]

	@pytest.mark.parametrize('to_kspec', [
		KmerSpec(10, 'ATGAC'),   # Reduce k
		KmerSpec(8, 'ATGAC'),    # Reduce k
		KmerSpec(9, 'ATGACGT'),  # Extend prefix
		KmerSpec(7, 'ATGACGT'),  # Extend prefix and reduce k further
	])
	def test_convert(self, seqs, to_kspec):
		from gambit.sigs.calc import calc_signature

		from_kspec = KmerSpec(11, 'ATGAC')

		for seq in seqs:
			from_sig = calc_signature(from_kspec, seq)
			from_vec = sparse_to_dense(from_kspec.k, from_sig)

			to_vec = convert_dense(from_kspec, to_kspec, from_vec)
			to_sig = convert_sparse(from_kspec, to_kspec, from_sig)

			found_sig = calc_signature(to_kspec, seq)

			assert np.array_equal(to_sig, found_sig)
			assert np.array_equal(to_vec, sparse_to_dense(to_kspec.k, found_sig))
