"""Store k-mer signature sets in HDF5 format."""

import json
from typing import Optional

import numpy as np
import h5py as h5

from .base import SignatureArray, ConcatenatedSignatureArray, AbstractSignatureArray, SignaturesMeta,\
	ReferenceSignatures
from gambit.kmers import KmerSpec
from gambit._cython.metric import BOUNDS_DTYPE
from gambit.util.io import FilePath


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

	Inherits from :class:`gambit.sigs.base.AbstractSignatureArray`, so behaves as a sequence of
	k-mer signatures supporting Numpy-style advanced indexing.

	Behaves as a context manager which yields itself on enter and closes the underlying HDF5 file
	object on exit. The :meth:`__bool__` method can be used to check whether the file is currently
	open and valid.

	Attributes
	----------
	group
		HDF5 group object data is read from.
	format_version
		Version of file format

	Parameters
	----------
	group
		Open, readable :class:`h5py.Group` or :class:`h5py.File` object.
	"""
	group: h5.Group
	format_version: int
	ids: h5.Dataset

	def __init__(self, group: h5.Group):
		self.group = group

		if FMT_VERSION_ATTR not in group.attrs:
			raise RuntimeError('HDF5 group does not contain a signature set')

		self.format_version = group.attrs[FMT_VERSION_ATTR]
		if self.format_version != CURRENT_FMT_VERSION:
			raise ValueError(f'Unrecognized format version: {self.format_version}')

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

	def close(self):
		"""Close the underlying HDF5 file."""
		if self.group:
			self.group.file.close()

	def __bool__(self):
		"""Check whether the underlying HDF5 file object is open."""
		return bool(self.group)

	def __enter__(self):
		return self

	def __exit__(self, *args):
		self.close()

	@classmethod
	def _init_attrs(cls, group: h5.Group, kmerspec: KmerSpec, meta: SignaturesMeta):
		"""Initialize attributes of group."""

		group.attrs[FMT_VERSION_ATTR] = CURRENT_FMT_VERSION
		group.attrs['kmerspec_k'] = kmerspec.k
		group.attrs['kmerspec_prefix'] = kmerspec.prefix_str

		write_metadata(group, meta)

	@classmethod
	def _init_datasets(cls, group: h5.Group, signatures: AbstractSignatureArray, ids: np.ndarray, values_kw = None):
		"""Initialize datasets of group."""

		if values_kw is None:
			values_kw = dict()

		if ids.dtype.kind == 'U':
			# h5py doesn't support writing Numpy U data type
			ids = ids.astype(object)
			ids_dtype = h5.string_dtype()
		elif ids.dtype.kind in 'OS':
			ids_dtype = h5.string_dtype()
		elif ids.dtype.kind in 'ui':
			ids_dtype = ids.dtype
		else:
			raise ValueError('ids array must contain integers or strings.')

		group.create_dataset('ids', data=ids, dtype=ids_dtype)

		if isinstance(signatures, SignatureArray):
			group.create_dataset('values', data=signatures.values, **values_kw)
			group.create_dataset('bounds', data=signatures.bounds, dtype=BOUNDS_DTYPE)

		else:
			n = len(signatures)
			sizes = np.asarray(signatures.sizes())

			bounds = group.create_dataset('bounds', shape=n + 1, dtype=BOUNDS_DTYPE)
			bounds[0] = 0
			bounds[1:] = np.cumsum(sizes, dtype=BOUNDS_DTYPE)

			values = group.create_dataset('values', shape=int(bounds[-1]), dtype=signatures.dtype, **values_kw)
			for i in range(n):
				values[bounds[i]:bounds[i + 1]] = signatures[i]

	@classmethod
	def create(cls,
	           group: h5.Group,
	           signatures: AbstractSignatureArray,
	           *,
	           compression: Optional[str] = None,
	           compression_opts = None,
	           ) -> 'HDF5Signatures':
		"""Store k-mer signatures and associated metadata in an HDF5 group.

		Parameters
		----------
		group
			HDF5 group to store data in.
		signatures
			Array of signatures to store. If an instance of
			:class:`gambit.sigs.base.ReferenceSignatures` its metadata will be stored as well,
			otherwise default/empty values will be used.
		compression
			Compression type for values array. One of ``['gzip', 'lzf', 'szip']``. See the
			section on
			`compression filters <https://docs.h5py.org/en/stable/high/dataset.html#lossless-compression-filters>`_
			in ``h5py``'s documentation.
		compression_opts
			Sets compression level (0-9) for gzip compression, no effect for other types.
		"""

		if isinstance(signatures, ReferenceSignatures):
			ids = np.asarray(signatures.ids)
			if ids.shape != (len(signatures),):
				raise ValueError('Length of ids must match length of data')

			meta = signatures.meta

		else:
			ids = np.arange(len(signatures))
			meta = SignaturesMeta()

		kw = dict(compression=compression, compression_opts=compression_opts)

		cls._init_attrs(group, signatures.kmerspec, meta)
		cls._init_datasets(group, signatures, ids, values_kw=kw)

		return cls(group)


def load_signatures_hdf5(path: FilePath, **kw) -> HDF5Signatures:
	"""Open HDF5 signature file.

	Parameters
	----------
	path
		File to open.
	\\**kw
		Additional keyword arguments to :func:`h5py.File`.
	"""
	return HDF5Signatures(h5.File(path, **kw))


def dump_signatures_hdf5(path: FilePath,
                         signatures: AbstractSignatureArray,
                         **kw,
                         ):
	"""Write k-mer signatures and associated metadata to an HDF5 file.

	Parameters
	----------
	path
		File to write to.
	signatures
		Array of signatures to store.
	\\**kw
		Additional keyword arguments to :meth:`HDF5Signatures.create`.
	"""
	with h5.File(path, 'w') as f:
		HDF5Signatures.create(f, signatures, **kw)
