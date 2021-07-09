"""Store k-mer signature sets in HDF5 format."""

import json
from typing import Sequence, Union, Optional

import numpy as np
import h5py as h5

from .base import SignaturesMeta, ReferenceSignatures
from .array import SignatureArray, ConcatenatedSignatureArray
from gambit.kmers import KmerSpec
from gambit.io.util import FilePath


#: Name of HDF5 group attribute which both stores the format version and also identifies the group
#: as containing signature data.
FMT_VERSION_ATTR = 'gambit_signatures_version'

#: Current version of the data format. Integer which should be incremented each time the format
#: changes.
CURRENT_FMT_VERSION = 1

STR_DTYPE = h5.string_dtype()


def none_to_empty(value, dtype: np.dtype):
	"""Convert None values to :class:`h5py.Empty`, passing other types through.
	"""
	return h5.Empty(dtype) if value is None else value

def empty_to_none(value):
	"""Convert :class:`h5py.Empty` instances to None, passing other types through.
	"""
	return None if isinstance(value, h5.Empty) else value


def write_metadata(group: h5.Group, meta: SignaturesMeta):
	"""Write signature set metadata to HDF5 group attributes."""
	group.attrs['id'] = none_to_empty(meta.id, STR_DTYPE)
	group.attrs['name'] = none_to_empty(meta.name, STR_DTYPE)
	group.attrs['id_attr'] = none_to_empty(meta.id_attr, STR_DTYPE)
	group.attrs['version'] = none_to_empty(meta.version, STR_DTYPE)
	group.attrs['description'] = none_to_empty(meta.description, STR_DTYPE)

	if meta.extra is not None:
		group.attrs['extra'] = json.dumps(meta.extra)
	else:
		group.attrs['extra'] = h5.Empty(STR_DTYPE)

def read_metadata(group: h5.Group) -> SignaturesMeta:
	"""Read signature set metadata from HDF5 group attributes."""
	extra_str = empty_to_none(group.attrs.get('extra'))
	extra = None if extra_str is None else json.loads(extra_str)

	return SignaturesMeta(
		id=empty_to_none(group.attrs.get('id')),
		name=empty_to_none(group.attrs.get('name')),
		id_attr=empty_to_none(group.attrs.get('id_attr')),
		version=empty_to_none(group.attrs.get('version')),
		description=empty_to_none(group.attrs.get('description')),
		extra=extra,
	)


class HDF5Signatures(ConcatenatedSignatureArray, ReferenceSignatures):
	"""Stores a set of k-mer signatures and associated metadata in an HDF5 group.

	Inherits from :class:`gambit.signatures.base.AbstractSignatureArray`, so behaves as a sequence of
	k-mer signatures supporting Numpy-style advanced indexing.

	Attributes
	----------
	group
		HDF5 group object data is read from.

	Parameters
	----------
	group
		Open, readable :class:`h5py.Group` or :class:`h5py.File` object.
	"""
	group: h5.Group
	ids: h5.Dataset

	def __init__(self, group: h5.Group):
		self.group = group

		if FMT_VERSION_ATTR not in group.attrs:
			raise RuntimeError('HDF5 group does not contain a signature set')

		version = group.attrs[FMT_VERSION_ATTR]
		if version != CURRENT_FMT_VERSION:
			raise ValueError(f'Unrecognized format version: {version}')

		self.kmerspec = KmerSpec(group.attrs['kmerspec_k'], group.attrs['kmerspec_prefix'])
		self.meta = read_metadata(group)

		self.values = group['values']
		self.bounds = group['bounds']

		ids_data = group['ids']
		if ids_data.dtype.kind == 'O':
			# String data set reads out bytes as default
			self.ids = ids_data.asstr()[:]
		else:
			self.ids = ids_data[:]

	@classmethod
	def _init_attrs(cls, group: h5.Group, kmerspec: KmerSpec, meta: SignaturesMeta):
		"""Initialize attributes of group."""

		group.attrs[FMT_VERSION_ATTR] = CURRENT_FMT_VERSION
		group.attrs['kmerspec_k'] = kmerspec.k
		group.attrs['kmerspec_prefix'] = kmerspec.prefix.decode('ascii')

		write_metadata(group, meta)

	@classmethod
	def _init_datasets(cls, group: h5.Group, data: SignatureArray, ids: np.ndarray):
		"""Initialize datasets of group."""

		if ids.dtype.kind == 'U':
			# h5py doesn't support writing Numpy U data type
			ids = ids.astype(object)
			ids_dtype = h5.string_dtype()
		elif ids.dtype.kind in 'OS':
			ids_dtype = h5.string_dtype()
		elif ids.dtype.kind in 'ui':
			ids_dtype=ids.dtype
		else:
			raise ValueError('ids array must contain integers or strings.')

		group.create_dataset('values', data=data.values)
		group.create_dataset('bounds', data=data.bounds)
		group.create_dataset('ids', data=ids, dtype=ids_dtype)

	@classmethod
	def open(cls, path: FilePath, **kw) -> 'HDF5Signatures':
		"""Open from file.

		Parameters
		----------
		path
			File to open.
		\\**kw
			Additional keyword arguments to :func:`h5py.File`.
		"""
		return cls(h5.File(path, **kw))

	@classmethod
	def create(cls,
	           group: h5.Group,
	           kmerspec: KmerSpec,
	           data: SignatureArray,
	           ids: Union[Sequence[int], Sequence[str], None] = None,
	           meta: Optional[SignaturesMeta] = None,
	           ) -> 'HDF5Signatures':
		"""Store k-mer signatures and associated metadata in an HDF5 group.

		Parameters
		----------
		group
			HDF5 group to store data in.
		kmerspec
			``KmerSpec`` used to calculate the signatures.
		data
			Array of signatures to store.
		ids
			Array of unique string or integer IDs for signatures in ``data``.  Defaults to
			consecutive integers starting from zero.
		meta
			Additional optional metadata to attach.
		"""

		if ids is None:
			ids = np.arange(len(data))
		else:
			ids = np.asarray(ids)
			if ids.shape != (len(data),):
				raise ValueError('Length of ids must match length of data')

		if meta is None:
			meta = SignaturesMeta()

		cls._init_attrs(group, kmerspec, meta)
		cls._init_datasets(group, data, ids)

		return cls(group)
