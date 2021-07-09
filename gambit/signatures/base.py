from typing import Optional, Sequence, Mapping, Any, Union
from abc import abstractmethod

import numpy as np
from attr import attrs, attrib

from gambit.kmers import KmerSpec, KmerSignature


@attrs()
class SignaturesMeta:
	"""Metadata describing a set of k-mer signatures.

	All attributes are optional.

	Attributes
	----------
	id
		Any kind of string ID that can be used to uniquely identify the signature set.
	version
		Version string (ideally PEP 440-compliant).
	name
		Short human-readable name.
	id_attr
		Name of ``Genome`` attribute the IDs correspond to (see :data:`gambit.db.models.GENOME_ID_ATTRS`).
		Optional, but signature set cannot be used as a reference for queries without it.
	description
		Human-readable description.
	extra
		Extra arbitrary metadata. Should be a ``dict`` or other mapping which can be converted to JSON.
	"""

	id : Optional[str] = attrib(default=None, kw_only=True)
	name : Optional[str] = attrib(default=None, kw_only=True)
	version : Optional[str] = attrib(default=None, kw_only=True)
	id_attr : Optional[str] = attrib(default=None, kw_only=True)
	description : Optional[str] = attrib(default=None, kw_only=True, repr=False)
	extra : Mapping[str, Any] = attrib(factory=dict, kw_only=True, repr=False)


def sigarray_eq(a1: Sequence, a2: Sequence) -> bool:
	"""Check two sequences of sparse k-mer signatures for equality."""
	return len(a1) == len(a2) and all(map(np.array_equal, a1, a2))


class AbstractSignatureArray(Sequence[KmerSignature]):
	"""
	Abstract base class for types which behave as a (non-mutable) sequence of k-mer signatures
	(k-mer sets in sparse coordinate format).

	Elements should be Numpy arrays with integer data type. Should implement numpy-style advanced
	indexing, see :class:`gambit.util.indexing.AdvancedIndexingMixin`. Slicing and advanced indexing
	should return another instance of ``AbstractSignatureArray``.

	Attributes
	----------
	dtype
		Numpy data type of signatures.
	"""
	dtype: np.dtype

	@abstractmethod
	def sizeof(self, index: int) -> int:
		"""Get the size/length of the signature at the given index.

		Should be the case that

		    sigarray.size_of(i) == len(sigarray[i])

		Parameters
		----------
		index
			Index of signature in array.
		"""

	def sizes(self) -> Sequence[int]:
		"""Get the sizes of all signatures in the array."""
		return np.fromiter(map(self.sizeof, range(len(self))))

	@abstractmethod
	def __getitem__(self, index: Union[int, slice, Sequence[int], Sequence[bool]]) -> Union[KmerSignature, 'AbstractSignatureArray']:
		pass

	def __eq__(self, other):
		if isinstance(other, Sequence):
			return sigarray_eq(self, other)
		else:
			return NotImplemented


class ReferenceSignatures(AbstractSignatureArray):
	"""Base class for an array of reference genome signatures plus metadata.

	This contains the extra data needed for the signatures to be used for running queries.

	Attributes
	----------
	kmerspec
		K-mer spec used to calculate signatures.
	ids
		Array of unique string or integer IDs for each signature. Length should be equal to length of
		``ReferenceSignatures`` object.
	meta
		Other metadata describing signatures.
	"""
	kmerspec: KmerSpec
	ids: Sequence
	meta: SignaturesMeta
