from typing import Sequence, Optional

import numpy as np

from .base import AbstractSignatureArray
from gambit._cython.metric import BOUNDS_DTYPE
from gambit.kmers import KmerSignature
from gambit.util.indexing import AdvancedIndexingMixin


class ConcatenatedSignatureArray(AdvancedIndexingMixin, AbstractSignatureArray):
	"""Base class for signature arrays which store signatures in a single data array.

	Attributes
	----------
	values
		K-mer signatures concatenated into single numpy-like array.
	bounds
		Numpy-like array storing indices bounding each individual k-mer signature in ``values``.
		The ``i``th signature is at ``values[bounds[i]:bounds[i + 1]]``.
	"""

	def __len__(self):
		return len(self.bounds) - 1

	def _getitem_int(self, i):
		return self.values[self.bounds[i]:self.bounds[i + 1]]

	def _getitem_slice(self, s):
		start, stop, step = s.indices(len(self))
		if step != 1 or stop <= start:
			return super()._getitem_slice(s)

		values = self.values[self.bounds[start]:self.bounds[stop]]
		bounds = self.bounds[start:(stop + 1)] - self.bounds[start]
		return SignatureArray.from_arrays(values, bounds)

	def _getitem_int_array(self, indices):
		out = SignatureArray.uninitialized([self.sizeof(i) for i in indices], dtype=self.values.dtype)
		for i, idx in enumerate(indices):
			np.copyto(out[i], self._getitem_int(idx), casting='unsafe')

		return out

	@property
	def dtype(self):
		return self.values.dtype

	def sizeof(self, index):
		i = self._check_index(index)
		return self.bounds[i + 1] - self.bounds[i]

	def sizes(self):
		return np.diff(self.bounds)


class SignatureArray(ConcatenatedSignatureArray):
	"""Stores a collection of k-mer signatures in a single contiguous Numpy array.

	This format enables the calculation of many Jaccard scores in parallel, see
	:func:`gambit.metric.jaccard_sparse_array`.

	Numpy-style indexing with an array of integers or bools is supported and will return another
	``SignatureArray``. If indexed with a contiguous slice the :attr:`values` of the returned
	array will be a view of the original instead of a copy.

	Attributes
	----------
	values
		K-mer signatures concatenated into single Numpy array.
	bounds
		Array storing indices bounding each individual k-mer signature in :attr:`values`.
		The ``i``th signature is at ``values[bounds[i]:bounds[i + 1]]``.
	"""
	values : np.ndarray
	bounds : np.ndarray

	@classmethod
	def _unint_arrays(cls, lengths, dtype):
		"""Get uninitialized values array and bounds array from signature lengths."""
		bounds = np.zeros(len(lengths) + 1, dtype=BOUNDS_DTYPE)
		np.cumsum(lengths, dtype=BOUNDS_DTYPE, out=bounds[1:])
		values = np.empty(bounds[-1], dtype=dtype)
		return values, bounds

	def _init_from_arrays(self, values, bounds):
		self.values = values
		self.bounds = bounds

	def __init__(self, signatures : Sequence[KmerSignature], dtype: Optional[np.dtype] = None):
		"""
		Parameters
		----------
		signatures
			Sequence of k-mer signatures.
		dtype
			Numpy dtype of :attr:`values` array.
		"""
		if isinstance(signatures, SignatureArray):
			# Can just copy arrays directly
			if dtype is None:
				values = signatures.values.copy()
			else:
				values = signatures.values.astype(dtype)
			bounds = signatures.bounds.copy()

			self._init_from_arrays(values, bounds)

		else:
			# Prepare with uninitialized values array
			lengths = list(map(len, signatures))
			values, bounds = self._unint_arrays(lengths, np.dtype('u8') if dtype is None else dtype)
			self._init_from_arrays(values, bounds)

			# Copy signatures to values array
			for i, sig in enumerate(signatures):
				np.copyto(self[i], sig, casting='unsafe')

	@classmethod
	def from_arrays(cls, values : np.ndarray, bounds : np.ndarray) -> 'SignatureArray':
		"""Create directly from values and bounds arrays."""
		sa = cls.__new__(cls)
		sa._init_from_arrays(values, bounds)
		return sa

	@classmethod
	def uninitialized(cls, lengths : Sequence[int], dtype: np.dtype = np.dtype('u8')) -> 'SignatureArray':
		"""Create with an uninitialized values array.

		Parameters
		----------
		lengths
			Sequence of lengths for each sub-array/signature.
		dtype
			Numpy dtype of shared coordinates array.
		"""
		return cls.from_arrays(*cls._unint_arrays(lengths, dtype))

	def __repr__(self):
		return f'<{type(self).__name__} length={len(self)} values.dtype={self.values.dtype}>'
