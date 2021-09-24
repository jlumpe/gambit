from abc import abstractmethod
from typing import NewType, Sequence, Union

import numpy as np

from gambit.kmers import KmerSpec

#: Type for k-mer signatures (k-mer sets in sparse coordinate format)
KmerSignature = NewType('KmerSignature', np.ndarray)
# TODO - use nptyping package to specify dimensions and data type?


def sigarray_eq(a1: Sequence, a2: Sequence) -> bool:
	"""Check two sequences of sparse k-mer signatures for equality."""
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
		return np.fromiter(map(self.sizeof, range(len(self))))

	@abstractmethod
	def __getitem__(self, index: Union[int, slice, Sequence[int], Sequence[bool]]) -> Union[KmerSignature, 'AbstractSignatureArray']:
		pass

	def __eq__(self, other):
		if isinstance(other, AbstractSignatureArray):
			return self.kmerspec == other.kmerspec and sigarray_eq(self, other)
		else:
			return NotImplemented
