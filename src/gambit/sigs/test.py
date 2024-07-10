"""Utilities for testing signature set types."""

import pytest
import numpy as np

from .base import AbstractSignatureArray, sigarray_eq


class AbstractSignatureArrayTests:
	"""Base for test classes which test :class:`AbstractSignatureArray` implementations.

	Must implement the following fixtures (unless they exist at the module level):

	* ``instance`` - ``AbstractSignatureArray`` instance to be tested.
	* ``ref_instance`` - Instance of reference implementation to compare to. May be a Numpy array of
	  object dtype containing the same signatures, or another type which is trusted to implement
	  ``AbstractSignatureArray`` indexing correctly.

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

	def check_getindex_scalar(self, instance, ref_instance, index, result, ref_result):
		"""Check indexing with a single integer.."""
		assert isinstance(result, np.ndarray)
		assert result.dtype == instance.dtype
		assert np.array_equal(result, ref_result)

	def check_getindex_subseq(self, instance, ref_instance, index, result, ref_result):
		"""Check any indexing behavior which results in a subsequence."""
		assert isinstance(result, AbstractSignatureArray)
		assert sigarray_eq(result, ref_result)
		assert result.kmerspec == instance.kmerspec
		assert result.dtype == instance.dtype

	def check_getindex_slice(self, instance, ref_instance, index, result, ref_result):
		"""Check indexing with a slice object.

		Defaults to calling :meth:`check_getindex_subseq`.
		"""
		return self.check_getindex_subseq(instance, ref_instance, index, result, ref_result)

	def check_getindex_int_array(self, instance, ref_instance, index, result, ref_result):
		"""Check indexing with an integer array/sequence.

		Defaults to calling :meth:`check_getindex_subseq`.
		"""
		return self.check_getindex_subseq(instance, ref_instance, index, result, ref_result)

	def check_getindex_bool_array(self, instance, ref_instance, index, result, ref_result):
		"""Check indexing with a boolean array/sequence.

		Defaults to calling :meth:`check_getindex_subseq`.
		"""
		return self.check_getindex_subseq(instance, ref_instance, index, result, ref_result)

	#
	# Fixtures
	#

	@pytest.fixture()
	def slice_indices(self, instance):
		"""Slice indices to test."""
		n = len(instance)
		return [
			slice(None, n // 2),
			slice(n // 3, 2 * n // 3),
			slice(None),
			slice(None, None, 1),
			slice(None, None, 2),
		]

	@pytest.fixture()
	def int_array_indices(self, instance):
		"""Integer array indices to test."""
		np.random.seed(0)

		indices = [
			[],
		]

		n = len(instance)
		if n > 0:
			# Some random indices
			index = np.random.choice(n, 10)

			# Add negative indices
			index = np.concatenate([index, index - n])

			indices.append(index)

		return indices

	@pytest.fixture()
	def bool_array_indices(self, instance):
		"""Boolean array indices to test."""
		return [
			np.arange(len(instance)) % 2 == 0,
		]

	#
	# Tests
	#

	def test_basic(self, instance, ref_instance):
		"""Test basic behavior and attributes."""

		# Check len
		assert len(instance) == len(ref_instance)

		# Check dtype attribute is np.dtype
		assert isinstance(instance.dtype, np.dtype)

	def test_sizes(self, instance, ref_instance):
		"""Test sizes() and sizeof() methods."""
		for i, sig in enumerate(ref_instance):
			assert instance.sizeof(i) == len(sig)
			assert instance.sizeof(np.int64(i)) == len(sig)

		with pytest.raises(IndexError):
			instance.sizeof(len(instance))

		assert np.array_equal(instance.sizes(), [instance.sizeof(i) for i in range(len(instance))])

	def test_iteration(self, instance, ref_instance):
		"""Test iteration protocol."""
		l = list(iter(instance))
		sigarray_eq(l, ref_instance)

	def test_getitem_single(self, instance, ref_instance):
		"""Test __getitem__ with a single integer."""
		n = len(ref_instance)

		for i in range(n):
			self.check_getindex_scalar(instance, ref_instance, i, instance[i], ref_instance[i])

	def test_getitem_slice(self, instance, ref_instance, slice_indices):
		"""Test __getitem__ with a slice object."""
		for s in slice_indices:
			self.check_getindex_slice(instance, ref_instance, s, instance[s], ref_instance[s])

	def test_getitem_int_array(self, instance, ref_instance, int_array_indices):
		"""Test __getitem__ with an integer array."""
		for index in int_array_indices:
			self.check_getindex_int_array(instance, ref_instance, index, instance[index], ref_instance[index])

	def test_getitem_bool_array(self, instance, ref_instance, bool_array_indices):
		"""Test __getitem__ with a boolean array."""
		for index in bool_array_indices:
			self.check_getindex_int_array(instance, ref_instance, index, instance[index], ref_instance[index])
