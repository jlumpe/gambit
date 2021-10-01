"""Array-like objects for storing k-mer signatures."""
from abc import abstractmethod
from typing import Sequence, Optional, Iterable, MutableSequence, Union

import numpy as np

from .base import KmerSignature
from gambit.kmers import KmerSpec
from gambit._cython.metric import BOUNDS_DTYPE
from gambit.util.indexing import AdvancedIndexingMixin


def sigarray_eq(a1: Sequence[KmerSignature], a2: Sequence[KmerSignature]) -> bool:
	"""Check two sequences of sparse k-mer signatures for equality.

	Unlike :meth:`.AbstractSignatureArray.__eq__` this works on any sequence type containing
	signatures and does not use the :attr:`.AbstractSignatureArray.kmerspec` attribute.
	"""
	return len(a1) == len(a2) and all(map(np.array_equal, a1, a2))


class AbstractSignatureArray(Sequence[KmerSignature]):
	"""
	Abstract base class for types which behave as a (non-mutable) sequence of k-mer signatures
	(k-mer sets in sparse coordinate format).

	The signature data itself may already be present in memory or may be loaded lazily from the file
	system when the object is indexed.

	Elements should be Numpy arrays with integer data type. Should implement numpy-style advanced
	indexing, see :class:`gambit.util.indexing.AdvancedIndexingMixin`. Slicing and advanced indexing
	should return another instance of ``AbstractSignatureArray``.

	Attributes
	----------
	kmerspec
		K-mer spec used to calculate signatures.
	dtype
		Numpy data type of signatures.
	"""
	kmerspec: KmerSpec
	dtype: np.dtype

	def sizeof(self, index: int) -> int:
		"""Get the size/length of the signature at the given index.

		Should be the case that

		    sigarray.size_of(i) == len(sigarray[i])

		Parameters
		----------
		index
			Index of signature in array.
		"""
		return len(self[index])

	def sizes(self) -> Sequence[int]:
		"""Get the sizes of all signatures in the array."""
		return np.fromiter(map(self.sizeof, range(len(self))), dtype=int)

	@abstractmethod
	def __getitem__(self, index: Union[int, slice, Sequence[int], Sequence[bool]]) -> Union[KmerSignature, 'AbstractSignatureArray']:
		pass

	def __eq__(self, other):
		"""Compare two ``AbstractSignatureArray`` instances for equality.

		Two instances are considered equal if they are equivalent as sequences (see
		:func:`.sigarray_eq`) and have the same :attr:`kmerspec`.
		"""
		if isinstance(other, AbstractSignatureArray):
			return self.kmerspec == other.kmerspec and sigarray_eq(self, other)
		else:
			return NotImplemented


class ConcatenatedSignatureArray(AdvancedIndexingMixin, AbstractSignatureArray):
	"""Base class for signature arrays which store signatures in a single data array.

	Attributes
	----------
	values
		K-mer signatures concatenated into single numpy-like array.
	bounds
		Numpy-like array storing indices bounding each individual k-mer signature in ``values``.
		The ``i``\\ th signature is at ``values[bounds[i]:bounds[i + 1]]``.
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
		return SignatureArray.from_arrays(values, bounds, self.kmerspec)

	def _getitem_int_array(self, indices):
		out = SignatureArray.uninitialized([self.sizeof(i) for i in indices], self.kmerspec, dtype=self.values.dtype)
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
	:func:`gambit.metric.jaccarddist_array`.

	Numpy-style indexing with an array of integers or bools is supported and will return another
	``SignatureArray``. If indexed with a contiguous slice the :attr:`values` of the returned
	array will be a view of the original instead of a copy.

	Attributes
	----------
	values
		K-mer signatures concatenated into single Numpy array.
	bounds
		Array storing indices bounding each individual k-mer signature in :attr:`values`.
		The ``i``\\ th signature is at ``values[bounds[i]:bounds[i + 1]]``.
	"""
	values : np.ndarray
	bounds : np.ndarray

	@classmethod
	def _uninit_arrays(cls, lengths, dtype):
		"""Get uninitialized values array and bounds array from signature lengths."""
		bounds = np.zeros(len(lengths) + 1, dtype=BOUNDS_DTYPE)
		np.cumsum(lengths, dtype=BOUNDS_DTYPE, out=bounds[1:])
		values = np.empty(bounds[-1], dtype=dtype)
		return values, bounds

	def _init_from_arrays(self, values, bounds, kmerspec):
		self.values = values
		self.bounds = bounds
		self.kmerspec = kmerspec

	def __init__(self, signatures: Sequence[KmerSignature], kmerspec: Optional[KmerSpec] = None, dtype: Optional[np.dtype] = None):
		"""
		Parameters
		----------
		signatures
			Sequence of k-mer signatures.
		kmerspec
			K-mer spec used to calculate signatures. If None will take from ``signatures`` if it is
			an :class:`AbstractSignatureArray` instance.
		dtype
			Numpy dtype of :attr:`values` array. If None will use dtype of first element of
			``signatures``.
		"""
		if kmerspec is None:
			if isinstance(signatures, AbstractSignatureArray):
				kmerspec = signatures.kmerspec
			else:
				raise TypeError('kmerspec cannot be None if signatures is not an instance of AbstractSignatureArray')

		if isinstance(signatures, SignatureArray):
			# Can just copy arrays directly
			if dtype is None:
				values = signatures.values.copy()
			else:
				values = signatures.values.astype(dtype)
			bounds = signatures.bounds.copy()

			self._init_from_arrays(values, bounds, kmerspec)

		else:
			# Prepare with uninitialized values array
			if dtype is None:
				# Get dtype from first signature
				dtype = signatures[0].dtype if signatures else kmerspec.index_dtype

			lengths = list(map(len, signatures))
			values, bounds = self._uninit_arrays(lengths, dtype)
			self._init_from_arrays(values, bounds, kmerspec)

			# Copy signatures to values array
			for i, sig in enumerate(signatures):
				np.copyto(self[i], sig, casting='unsafe')

	@classmethod
	def from_arrays(cls, values: np.ndarray, bounds: np.ndarray, kmerspec: KmerSpec) -> 'SignatureArray':
		"""Create directly from values and bounds arrays."""
		sa = cls.__new__(cls)
		sa._init_from_arrays(values, bounds, kmerspec)
		return sa

	@classmethod
	def uninitialized(cls, lengths: Sequence[int], kmerspec: KmerSpec, dtype: np.dtype = None) -> 'SignatureArray':
		"""Create with an uninitialized values array.

		Parameters
		----------
		lengths
			Sequence of lengths for each sub-array/signature.
		kmerspec
		dtype
			Numpy dtype of shared coordinates array.
		"""
		values, bounds = cls._uninit_arrays(lengths, kmerspec.index_dtype if dtype is None else dtype)
		return cls.from_arrays(values, bounds, kmerspec)

	def __repr__(self):
		return f'<{type(self).__name__} length={len(self)} values.dtype={self.values.dtype}>'


class SignatureList(AdvancedIndexingMixin, AbstractSignatureArray, MutableSequence[KmerSignature]):
	"""Stores a collection of k-mer signatures in a standard Python list.

	Compared to :class:`SignatureArray` this isn't as efficient to calculate Jaccard scores with,
	but supports mutation and won't have to copy signatures to a new array on creation.
	"""

	def __init__(self, signatures: Iterable[KmerSignature], kmerspec: Optional[KmerSpec] = None, dtype: Optional[np.dtype] = None):
		"""
		Parameters
		----------
		signatures
			Iterable of k-mer signatures.
		kmerspec
			K-mer spec used to calculate signatures. If None will take from ``signatures`` if it is
			an :class:`AbstractSignatureArray` instance.
		dtype
			Numpy dtype of signatures. If None will use dtype of first element of
			``signatures``.
		"""
		self._list = list(signatures)

		if kmerspec is not None:
			self.kmerspec = kmerspec
		elif isinstance(signatures, AbstractSignatureArray):
			self.kmerspec = signatures.kmerspec
		else:
			raise TypeError('kmerspec cannot be None if signatures is not an instance of AbstractSignatureArray')

		if dtype is not None:
			self.dtype = dtype
		elif isinstance(signatures, AbstractSignatureArray):
			self.dtype = signatures.dtype
		elif len(self._list) > 0:
			self.dtype = self._list[0].dtype
		else:
			self.dtype = self.kmerspec.index_dtype

	def __len__(self):
		return len(self._list)

	def __iter__(self):
		return iter(self._list)

	def _getitem_int(self, i: int):
		return self._list[i]

	def _getitem_int_array(self, indices: np.ndarray):
		return SignatureList([self._list[i] for i in indices], self.kmerspec, self.dtype)

	def __setitem__(self, i: int, sig: KmerSignature):
		self._list[i] = sig

	def __delitem__(self, i: int):
		del self._list[i]

	def insert(self, i: int, sig: KmerSignature):
		self._list.insert(i, sig)
