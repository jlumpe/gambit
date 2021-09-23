"""Store metadata with k-mer signatures."""
from typing import Optional, Mapping, Any, Sequence

from attr import attrs, attrib

from gambit.kmers import KmerSpec
from .base import AbstractSignatureArray


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


class ReferenceSignatures(AbstractSignatureArray):
	"""Base class for an array of reference genome signatures plus metadata.

	This contains the extra data needed for the signatures to be used for running queries.

	Attributes
	----------
	ids
		Array of unique string or integer IDs for each signature. Length should be equal to length of
		``ReferenceSignatures`` object.
	meta
		Other metadata describing signatures.
	"""
	ids: Sequence
	meta: SignaturesMeta
