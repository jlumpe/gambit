"""Test gambit.signatures.SignatureArray."""

import pytest
import numpy as np

from gambit.signatures import SignatureArray, SignatureList
from gambit.kmers import KmerSpec
from gambit.test import make_signatures
from gambit.signatures.test import AbstractSignatureArrayTests


@pytest.fixture(params=[None, 'i8', 'u4'])
def sigarray(request):
	np.random.seed(0)
	return make_signatures(8, 100, request.param)


class TestSignatureArray:
	"""Test the SignatureArray class."""

	class TestAbstractSignatureArrayImplementation(AbstractSignatureArrayTests):
		"""Test implementation of the AbstractSignatureArray interface."""

		@pytest.fixture()
		def instance(self, sigarray):
			"""Instance to test."""
			return sigarray

		@pytest.fixture()
		def ref_instance(self, sigarray):
			"""Numpy array equivalent to `sigarray`."""
			return np.asarray(sigarray, dtype=object)

		def check_getindex_subseq(self, instance, ref_instance, index, result, ref_result):
			super().check_getindex_subseq(instance, ref_instance, index, result, ref_result)
			assert isinstance(result, SignatureArray)

	def test_uninitialized(self, sigarray):
		"""Test creating with the uninitialized() class method."""

		sa2 = SignatureArray.uninitialized(sigarray.sizes(), sigarray.kmerspec)
		assert len(sa2) == len(sigarray)
		assert np.array_equal(sa2.sizes(), sigarray.sizes())
		assert sa2.dtype == sigarray.kmerspec.index_dtype

		for i, sig in enumerate(sa2):
			assert sa2.sizeof(i) == len(sig)

	def test_construct_from_list(self, sigarray):
		"""Test construction from generic sequence type containing signatures."""

		sa2 = SignatureArray(list(sigarray), sigarray.kmerspec)
		assert sa2.dtype == sigarray.dtype
		assert sa2 == sigarray

		# Set different dtype
		for dtype in map(np.dtype, ['i8', 'u4']):
			sa2 = SignatureArray(list(sigarray), sigarray.kmerspec, dtype=dtype)
			assert sa2.dtype == dtype
			assert sa2 == sigarray

	def test_construct_from_signaturearray(self, sigarray):
		"""Test construction from another SignatureArray."""

		sa2 = SignatureArray(sigarray)
		assert sa2.kmerspec == sigarray.kmerspec
		assert sa2.dtype == sigarray.dtype
		assert sa2 == sigarray

		# Set different dtype
		for dtype in map(np.dtype, ['i8', 'u4']):
			sa2 = SignatureArray(sigarray, dtype=dtype)
			assert sa2.dtype == dtype
			assert sa2 == sigarray

		# Set different kmerspec
		kspec = KmerSpec(sigarray.kmerspec.k, sigarray.kmerspec.prefix_str + 'A')
		sa2 = SignatureArray(sigarray, kspec)
		assert sa2.kmerspec == kspec
		assert sa2 != sigarray

	def test_empty(self):
		"""Really an edge case, but test it anyways."""

		sigarray = SignatureArray([], KmerSpec(1, 'A'))
		assert len(sigarray) == 0

		with pytest.raises(IndexError):
			sigarray[0]


class TestSignatureList:
	"""Test the SignatureList class."""

	class TestAbstractSignatureArrayImplementation(AbstractSignatureArrayTests):
		"""Test implementation of AbstractSignatureArray interface."""

		@pytest.fixture()
		def instance(self, sigarray):
			"""Instance to test."""
			return SignatureList(sigarray)

		@pytest.fixture()
		def ref_instance(self, sigarray):
			"""Numpy array equivalent to `sigarray`."""
			return np.asarray(sigarray, dtype=object)

		def check_getindex_subseq(self, instance, ref_instance, index, result, ref_result):
			super().check_getindex_subseq(instance, ref_instance, index, result, ref_result)
			assert isinstance(result, SignatureList)

	def test_construct_from_list(self, sigarray):
		"""Test construction from generic sequence type containing signatures."""

		sl = SignatureList(list(sigarray), sigarray.kmerspec)
		assert sl.dtype == sigarray.dtype
		assert sl == sigarray

		# Set different dtype
		for dtype in map(np.dtype, ['i8', 'u4']):
			sl = SignatureList(list(sigarray), sigarray.kmerspec, dtype=dtype)
			assert sl.dtype == dtype
			assert sl == sigarray

	def test_construct_from_signaturearray(self, sigarray):
		"""Test construction from another AbstractSignatureArray instance."""

		sl = SignatureList(sigarray)
		assert sl.kmerspec == sigarray.kmerspec
		assert sl.dtype == sigarray.dtype
		assert sl == sigarray

		# Set different dtype
		for dtype in map(np.dtype, ['i8', 'u4']):
			sl = SignatureList(sigarray, dtype=dtype)
			assert sl.dtype == dtype
			assert sl == sigarray

		# Set different kmerspec
		kspec = KmerSpec(sigarray.kmerspec.k, sigarray.kmerspec.prefix_str + 'A')
		sl = SignatureList(sigarray, kspec)
		assert sl.kmerspec == kspec
		assert sl != sigarray

	def test_empty(self, sigarray):
		"""Test construction from empty list."""
		sl = SignatureList([], sigarray.kmerspec)
		assert sl.dtype == sigarray.kmerspec.index_dtype
