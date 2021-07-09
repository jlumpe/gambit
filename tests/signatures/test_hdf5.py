"""Test gambit.signatures.hdf5."""

import pytest
import h5py as h5
import numpy as np

from gambit.signatures.hdf5 import HDF5Signatures, read_metadata, write_metadata
from gambit.signatures import SignaturesMeta
from gambit.signatures.test import AbstractSignatureArrayTests
from gambit.kmers import KmerSpec
from gambit.test import make_signatures


# JSON data to use for metadata extra field
EXTRA = dict(
	num=4,
	string='foo',
	bool=True,
	array=[1, 2, 3],
)


@pytest.mark.parametrize('optional_attrs', [False, True])
# def test_metadata(h5file, optional_attrs):
def test_metadata(tmp_path, optional_attrs):
	"""Test reading/writing metadata"""

	fname = tmp_path / 'test.h5'

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


class TestHDF5Signatures:

	@pytest.fixture(scope='class')
	def kspec(self):
		return KmerSpec(8, 'ATG')

	@pytest.fixture(scope='class', params=[(1000, 'u8'), (1000, 'i4'), (0, 'u8')])
	def sigs(self, request, kspec):
		n, dtype = request.param
		return make_signatures(kspec.k, n, dtype)

	@pytest.fixture(scope='class')
	def meta(self):
		return SignaturesMeta(
			id='test',
			extra=EXTRA,
		)

	@pytest.fixture(scope='class', params=[int, str])
	def sig_ids(self, request, sigs):
		if request.param is int:
			return np.arange(len(sigs))
		elif request.param is str:
			return [f'test-{i}' for i in range(len(sigs))]
		else:
			assert 0

	@pytest.fixture(scope='class')
	def h5file(self, tmp_path_factory, sigs, kspec, sig_ids, meta):
		"""Write signatures to file and return file name."""
		fname = tmp_path_factory.mktemp('HDF5FileSignatures') / 'test.h5'

		with h5.File(fname, 'w') as f:
			HDF5Signatures.create(f, kspec, sigs, ids=sig_ids, meta=meta)

		return fname

	@pytest.fixture()
	def h5sigs(self, h5file):
		"""Open HDF5Signatures object."""
		with h5.File(h5file, 'r') as f:
			yield HDF5Signatures(f)

	def test_attrs(self, h5sigs, sigs, kspec, sig_ids, meta):
		"""Test basic attributes."""
		assert h5sigs.kmerspec == kspec
		assert h5sigs.dtype == sigs.values.dtype
		assert np.array_equal(h5sigs.ids, sig_ids)
		assert h5sigs.meta == meta

	class TestAbstractSignatureArrayImplementation(AbstractSignatureArrayTests):
		"""Test implementation of AbstractSignatureArray."""

		@pytest.fixture()
		def sigarray(self, h5sigs):
			return h5sigs

		@pytest.fixture()
		def refarray(self, sigs):
			return sigs
