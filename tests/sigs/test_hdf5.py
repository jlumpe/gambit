"""Test gambit.sigs.hdf5."""

from pathlib import Path

import pytest
import h5py as h5
import numpy as np

from gambit.sigs.hdf5 import read_metadata, write_metadata, load_signatures_hdf5, \
	dump_signatures_hdf5, HDF5Signatures
from gambit.sigs.base import SignaturesMeta, SignatureList, AnnotatedSignatures, \
	AbstractSignatureArray, SignaturesFileError, SignatureArray
from gambit.kmers import KmerSpec
from ..common import make_signatures
from .common import AbstractSignatureArrayTests


# JSON data to use for metadata extra field
EXTRA = dict(
	num=4,
	string='foo',
	bool=True,
	array=[1, 2, 3],
)


@pytest.mark.parametrize('optional_attrs', [False, True])
def test_metadata(tmp_path: Path, optional_attrs: bool):
	"""Test reading/writing metadata"""

	fname = tmp_path / 'test.gs'

	meta = SignaturesMeta(
		id='test',
		name='test',
		id_attr='refseq_acc' if optional_attrs else None,
		version='1.0' if optional_attrs else None,
		description='Test metadata' if optional_attrs else None,
		extra=EXTRA if optional_attrs else None,
	)

	with h5.File(fname, 'w') as file:
		write_metadata(file, meta)

	with h5.File(fname) as file:
		meta2 = read_metadata(file)

	assert meta2 == meta


def dump_load(sigs: AbstractSignatureArray, path: Path, **kw) -> HDF5Signatures:
	"""Dump signatures to HDF5 file and load them again."""
	f = path / 'test.gs'
	dump_signatures_hdf5(f, sigs, **kw)
	return load_signatures_hdf5(f)


def test_open_not_hdf5(tmp_path: Path):
	"""Test opening an invalid file."""

	# Not an HDF5 file
	file = tmp_path / 'not-hdf5.gs'
	with open(file, 'w') as f:
		f.write('foo')

	with pytest.raises(SignaturesFileError) as einfo:
		load_signatures_hdf5(file)

	assert einfo.value.filename == str(file)
	assert einfo.value.format == 'hdf5'


def test_open_invalid(tmp_path: Path):
	"""Test opening an invalid HDF5 file."""

	file = tmp_path / 'invalid.gs'
	with h5.File(file, 'w') as f:
		pass  # Empty

	with pytest.raises(SignaturesFileError) as einfo:
		load_signatures_hdf5(file)

	assert einfo.value.filename == str(file)
	assert einfo.value.format == 'hdf5'


class TestHDF5Signatures:

	@pytest.fixture(scope='class')
	def kspec(self):
		return KmerSpec(8, 'ATG')

	@pytest.fixture(scope='class', params=[(1000, 'u8'), (1000, 'i4'), (0, 'u8')])
	def sigs(self, request, kspec: KmerSpec):
		n, dtype = request.param
		return make_signatures(kspec, n, dtype)

	@pytest.fixture(scope='class')
	def h5file(self, tmp_path_factory, sigs: SignatureArray):
		"""Write signatures to file and return file name."""
		fname = tmp_path_factory.mktemp('HDF5FileSignatures') / 'test.gs'
		dump_signatures_hdf5(fname, sigs)
		return fname

	@pytest.fixture()
	def h5sigs(self, h5file: Path):
		"""Open HDF5Signatures object."""
		with load_signatures_hdf5(h5file) as sigs:
			yield sigs

	def test_attrs(self, h5sigs: AnnotatedSignatures, sigs: SignatureArray):
		"""Test basic attributes for signatures saved without metadata."""
		assert h5sigs.kmerspec == sigs.kmerspec
		assert h5sigs.dtype == sigs.values.dtype
		assert np.array_equal(h5sigs.ids, np.arange(len(h5sigs)))
		assert h5sigs.meta == SignaturesMeta()

	@pytest.mark.parametrize('id_type', [int, str])
	def test_attrs_meta(self, sigs: SignatureArray, id_type: type, tmp_path: Path):
		"""Test basic attributes for signatures saved with metadata."""

		if id_type is int:
			ids = np.arange(len(sigs)) + 1
		elif id_type is str:
			ids = [f'test-{i+1}' for i in range(len(sigs))]
		else:
			assert 0

		meta = SignaturesMeta(
			id='test',
			extra=EXTRA,
		)

		annotated = AnnotatedSignatures(sigs, ids, meta)

		with dump_load(annotated, tmp_path) as h5sigs:
			assert h5sigs.kmerspec == sigs.kmerspec
			assert h5sigs.dtype == sigs.dtype
			assert np.array_equal(h5sigs.ids, ids)
			assert h5sigs.meta == meta

	def test_close(self, h5sigs: HDF5Signatures):
		assert h5sigs.group
		assert h5sigs

		h5sigs.close()
		assert not h5sigs.group
		assert not h5sigs

		h5sigs.close()

	def test_context(self, h5sigs: HDF5Signatures):
		with h5sigs as value:
			assert value is h5sigs
			assert h5sigs.group
			assert h5sigs

		assert not h5sigs.group
		assert not h5sigs

	def test_create_from_list(self, sigs, tmp_path: Path):
		"""Test creating from other AbstractSignatureArray type."""
		siglist = SignatureList(sigs)

		with dump_load(siglist, tmp_path) as h5sigs:
			assert h5sigs == siglist

	@pytest.mark.parametrize('from_list', [False, True])
	@pytest.mark.parametrize('compression_level', [None, 7])
	def test_compression(self, from_list: bool, compression_level, sigs: SignatureArray, tmp_path: Path):
		"""Test creating with gzip compression."""
		create_from = SignatureList(sigs) if from_list else sigs

		with dump_load(create_from, tmp_path) as h5sigs:
			assert h5sigs == sigs

	class TestAbstractSignatureArrayImplementation(AbstractSignatureArrayTests):
		"""Test implementation of AbstractSignatureArray."""

		@pytest.fixture()
		def instance(self, h5sigs):
			return h5sigs

		@pytest.fixture()
		def ref_instance(self, sigs):
			return sigs
