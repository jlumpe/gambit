"""Utilities for testing signature set types."""

import pytest
import numpy as np

from .array import AbstractSignatureArray
from .base import sigarray_eq


class AbstractSignatureArrayTests:
	"""Base for test classes which test :class:`gambit.sigatures.base.AbstractSignatureArray` implementations.

	Must implement the following fixtures (unless they exist at the module level):
	* `sigarray` - ``AbstractSignatureArray`` instance to be tested.
	* `refarray` - Instance of reference implementation to compare to. May be a Numpy array of object
	  dtype containing the same signatures, or another type which is trusted to implement the
	  ``AbstractSignatureArray`` interface correctly.

	Methods which can be overridden:
	* :meth:`check_scalar_result`
	* :meth:`check_subseq_result`
	* :meth:`check_int_array_result`
	* :meth:`check_bool_array_result`
	* :meth:`check_slice_result`

	Fixtures which can be overridden:
	* ``slice_indices``
	* ``int_array_indices``
	* ``bool_array_indices``
	"""

	#
	# Helper methods
	#

	def check_getindex_scalar(self, sigarray, refarray, index, result, refresult):
		"""Check indexing with a single integer.

		By default just compares results based on :meth:`signature_eq`.
		"""
		assert np.array_equal(result, refresult)

	def check_getindex_subseq(self, sigarray, refarray, index, result, refresult):
		"""Check  any indexing behavior which results in a subsequence.

		By default just compares results based on :meth:`sigarray_eq`.
		"""
		assert isinstance(result, AbstractSignatureArray)
		assert sigarray_eq(result, refresult)
		assert result.dtype == sigarray.dtype

	def check_getindex_slice(self, sigarray, refarray, index, result, refresult):
		"""Check indexing with a slice object.

		Defaults to calling :meth:`check_getindex_subseq`.
		"""
		return self.check_getindex_subseq(sigarray, refarray, index, result, refresult)

	def check_getindex_int_array(self, sigarray, refarray, index, result, refresult):
		"""Check indexing with an integer array/sequence.

		Defaults to calling :meth:`check_getindex_subseq`.
		"""
		return self.check_getindex_subseq(sigarray, refarray, index, result, refresult)

	def check_getindex_bool_array(self, sigarray, refarray, index, result, refresult):
		"""Check indexing with a boolean array/sequence.

		Defaults to calling :meth:`check_getindex_subseq`.
		"""
		return self.check_getindex_subseq(sigarray, refarray, index, result, refresult)

	#
	# Fixtures
	#

	@pytest.fixture()
	def slice_indices(self, sigarray):
		"""Slice indices to test."""
		n = len(sigarray)
		return [
			slice(None, n // 2),
			slice(n // 3, 2 * n // 3),
			slice(None),
			slice(None, None, 1),
			slice(None, None, 2),
		]

	@pytest.fixture()
	def int_array_indices(self, sigarray):
		"""Integer array indices to test."""
		np.random.seed(0)

		indices = [
			[],
		]

		n = len(sigarray)
		if n > 0:
			# Some random indices
			index = np.random.choice(n, 10)

			# Add negative indices
			index = np.concatenate([index, index - n])

			indices.append(index)

		return indices

	@pytest.fixture()
	def bool_array_indices(self, sigarray):
		"""Boolean array indices to test."""
		return [
			np.arange(len(sigarray)) % 2 == 0,
		]

	#
	# Tests
	#

	def test_basic(self, sigarray, refarray):
		"""Test basic behavior and attributes."""

		# Check len
		assert len(sigarray) == len(refarray)

		# Check dtype attribute is np.dtype
		assert isinstance(sigarray.dtype, np.dtype)

	def test_sizes(self, sigarray, refarray):
		"""Test sizes() and sizeof() methods."""
		for i, sig in enumerate(refarray):
			assert sigarray.sizeof(i) == len(sig)
			assert sigarray.sizeof(np.int64(i)) == len(sig)

		with pytest.raises(IndexError):
			sigarray.sizeof(len(sigarray))

		assert np.array_equal(sigarray.sizes(), [sigarray.sizeof(i) for i in range(len(sigarray))])

	def test_iteration(self, sigarray, refarray):
		"""Test iteration protocol."""
		l = list(iter(sigarray))
		sigarray_eq(l, refarray)

	def test_getitem_single(self, sigarray, refarray):
		"""Test __getitem__ with a single integer."""
		n = len(refarray)

		for i in range(n):
			self.check_getindex_scalar(sigarray, refarray, i, sigarray[i], refarray[i])

	def test_getitem_slice(self, sigarray, refarray, slice_indices):
		"""Test __getitem__ with a slice object."""
		for s in slice_indices:
			self.check_getindex_slice(sigarray, refarray, s, sigarray[s], refarray[s])

	def test_getitem_int_array(self, sigarray, refarray, int_array_indices):
		"""Test __getitem__ with an integer array."""
		for index in int_array_indices:
			self.check_getindex_int_array(sigarray, refarray, index, sigarray[index], refarray[index])

	def test_getitem_bool_array(self, sigarray, refarray, bool_array_indices):
		"""Test __getitem__ with a boolean array."""
		for index in bool_array_indices:
			self.check_getindex_int_array(sigarray, refarray, index, sigarray[index], refarray[index])
