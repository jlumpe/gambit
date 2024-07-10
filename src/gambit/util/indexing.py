"""Tools to implement and test Numpy-style advanced indexing for Sequence types."""

from typing import Sequence

import numpy as np


class AdvancedIndexingMixin:
	"""Mixin class for sequence types which enables Numpy-style advanced indexing.

	Implements a ``__getitem__`` which takes care of identifying the type of index, performing
	bounds checking, and converting negative indices.

	The following methods must be implemented by subtypes:
	* :meth:`_getitem_int`
	* :meth:`_getitem_int_array`

	The following methods may optionally be overridden, but default to calling :meth:`_getitem_int_array`:
	* :meth:`_getitem_range`
	* :meth:`_getitem_bool_array`
	"""

	def _check_index(self, i: int) -> int:
		"""Check integer index is in bounds, converting negative indices to positive."""
		i2 = i + len(self) if i < 0 else i
		if not 0 <= i2 < len(self):
			raise IndexError(f'Index out of bounds: {i}')
		return i2

	def _getitem_int(self, i: int):
		"""Get single value given valid integer index.

		Parameters
		----------
		i
			Positive in-bounds integer index.
		"""
		raise NotImplementedError()

	def _getitem_slice(self, index: slice):
		"""Get subsequence by indexing with a slice.

		Parameters
		----------
		slice
			Slice object. Guaranteed to be valid (start, stop, and step are integers or None and
			step size is non-zero).
		"""
		start, stop, step = index.indices(len(self))
		return self._getitem_int_array(np.arange(start, stop, step))

	def _getitem_int_array(self, index: Sequence[int]):
		"""Get subsequence given an array of integer indices.

		Parameters
		----------
		index
			Array of integer indices, all positive and in-bounds.
		"""
		raise NotImplementedError()

	def _getitem_bool_array(self, index: Sequence[bool]):
		"""Get subsequence given a boolean array.

		Parameters
		----------
		index
			Boolean Numpy array of same length as sequence.
		"""
		return self._getitem_int_array(np.flatnonzero(index))

	def __getitem__(self, index):
		input_index = index

		# Single integer index
		if isinstance(index, (int, np.integer)):
			return self._getitem_int(self._check_index(index))

		# Slice
		elif isinstance(index, slice):
			for i in [index.start, index.stop, index.step]:
				if i is not None and not isinstance(i, (int, np.integer)):
					raise TypeError('Slice indices must be integers or None')

			if index.step == 0:
				raise ValueError('Slice step cannot be zero')

			return self._getitem_slice(index)

		# Otherwise assume sequence of ints or bools, use Numpy to figure out array interpretation
		elif not isinstance(index, np.ndarray):
			# Special case - if an empty sequence np.asarray won't be able to infer dtype and
			# will default to floats
			if len(index) == 0:
				index = np.empty(0, dtype=int)

			else:
				try:
					index = np.asarray(index)
				except Exception as e:
					raise IndexError('Indices must be integers, slices, or integer or boolean sequences.') from e

		if index.ndim != 1:
			raise IndexError('Index arrays must be one-dimensional.')

		# Boolean array
		if index.dtype.kind == 'b':
			if len(index) != len(self):
				raise IndexError('Length of boolean index array does not match length of sequence.')
			return self._getitem_bool_array(index)

		# Integer array
		elif index.dtype.kind in 'iu':
			# Check bounds
			for i in index:
				self._check_index(i)

			# Convert negative indices to positive
			isneg = index < 0
			if isneg.any():
				# Don't modify input array
				if index is input_index:
					index = index.copy()
				np.add(index, len(self), out=index, where=isneg)

			return self._getitem_int_array(index)

		# Invalid dtype
		else:
			raise IndexError('Index arrays must have integer or boolean data type.')
