"""Test gambit.util.indexing module."""

from collections.abc import Sequence

import pytest
import numpy as np

from gambit.util.indexing import AdvancedIndexingMixin


class AdvancedIndexingTestClass(AdvancedIndexingMixin, Sequence):
	"""Just wraps a Numpy array."""

	def __init__(self, array):
		self.array = array
		self.calls = []  # Tracks which private functions were called in last call to __getitem__

	def __len__(self):
		return len(self.array)

	def _getitem_int(self, i):
		self.calls.append('int')
		assert isinstance(i, (int, np.integer))
		assert 0 <= i < len(self)
		return self.array[i]

	def _getitem_slice(self, index):
		self.calls.append('slice')
		assert slice.step != 0
		return super()._getitem_slice(index)

	def _getitem_int_array(self, index):
		self.calls.append('int_array')

		assert isinstance(index, np.ndarray)
		assert index.ndim == 1
		assert index.dtype.kind in 'iu'

		return self.array[index]

	def _getitem_bool_array(self, index):
		self.calls.append('bool_array')

		assert isinstance(index, np.ndarray)
		assert index.dtype.kind == 'b'
		assert index.shape == (len(self),)

		return super()._getitem_bool_array(index)

	def __getitem__(self, index):
		self.calls = []
		return super().__getitem__(index)


class TestAdvancedIndexingMixin:
	"""Test basic sequence class using AdvancedIndexingMixin."""

	@pytest.fixture(params=[False, True])
	def seq(self, request):
		if request.param:
			np.random.seed(0)
			arr = np.random.choice(1000, 100)
		else:
			# Make sure all tests work on empty sequence
			arr = np.empty(0, dtype=int)

		return AdvancedIndexingTestClass(arr)

	def test_iteration(self, seq):
		"""Test iterating over the sequence."""
		i = 0

		for x in seq:
			assert x == seq.array[i]
			i += 1

		assert i == len(seq)

	def test_getitem_single(self, seq):
		"""Test indexing with single integers."""
		n = len(seq)

		# All valid indices (positive and negative)
		indices = list(range(n))
		indices.extend(range(-1, -(n + 1), -1))

		for i in indices:
			for index in [i, np.int64(i)]:
				assert seq[index] == seq.array[i]
				assert seq.calls == ['int']

		# Out-of-bounds
		for i in [n, -(n + 1)]:
			for index in [i, np.int64(i)]:
				with pytest.raises(IndexError):
					seq[index]

	def test_getitem_slice(self, seq):
		"""Test indexing with slices."""
		n = len(seq)
		n_13 = n // 3
		n_23 = n_13 * 2

		slices = [
			slice(None),  # All
			slice(n_13),  # Start to index
			slice(n_23, None),  # Index to end
			slice(n_13, n_23),  # Index to index
			slice(n_13 - n, n_23 - n),  # Negative indices
			slice(n * 2),  # Past end, clipped to (length - 1)
			slice(-n * 2, None),  # Past beginning, clipped to 0
			slice(n_23, n_13),  # Empty
			slice(None, None, 1),  # Same as slice(None)
			slice(None, None, 2),  # Every other
			slice(n_13, None, -1),  # Backwards from index to start
			slice(None, n_13, 2),  # Every other with start index
		]

		for s in slices:
			assert np.array_equal(seq[s], seq.array[s])
			assert seq.calls == ['slice', 'int_array']

		# Zero step
		with pytest.raises(ValueError):
			seq[::0]

	def test_getitem_int_array(self, seq):
		"""Test indexing with integer arrays."""
		n = len(seq)

		indices = [
			[],
		]

		if n > 0:
			np.random.seed(0)
			index = np.random.choice(n, 10)  # 10 random in-bounds integers
			index = np.concatenate([index, index - n])  # Include negative versions
			indices.append(list(index))

		for indexlist in indices:
			for index in [indexlist, np.asarray(indexlist, dtype=int)]:
				assert np.array_equal(seq[index], seq.array[index])
				assert seq.calls == ['int_array']

		# Test with out-of-bounds indices
		for indexlist in indices:
			for i in [len(seq), -(len(seq) + 1)]:
				index = [*indexlist, i]
				with pytest.raises(IndexError):
					seq[index]

				index = np.asarray(index)
				with pytest.raises(IndexError):
					seq[index]

	def test_getitem_bool_array(self, seq):
		"""Test indexing with boolean arrays."""
		indices = [
			np.arange(len(seq)) % 2 == 0,
		]

		for index in indices:
			assert np.array_equal(seq[index], seq.array[index])
			assert seq.calls == ['bool_array', 'int_array']

			# Wrong length
			index2 = np.concatenate([index, [True]])
			with pytest.raises(IndexError):
				seq[index2]
