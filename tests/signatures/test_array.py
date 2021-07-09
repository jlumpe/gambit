"""Test gambit.signatures.SignatureArray."""

import pytest
import numpy as np

from gambit.signatures import SignatureArray
from gambit.test import make_signatures
from gambit.signatures.test import AbstractSignatureArrayTests


@pytest.fixture(params=[None, 'i8', 'u4'])
def sigarray(request):
	np.random.seed(0)
	return make_signatures(8, 100, request.param)


@pytest.fixture()
def refarray(sigarray):
	"""Numpy array equivalent to `sigarray`."""
	return np.asarray(sigarray, dtype=object)


class TestAbstractSignatureArrayImplementation(AbstractSignatureArrayTests):
	"""Test implementation of AbstractSignatureArray interface."""

	def check_getindex_scalar(self, sigarray, refarray, index, result, refresult):
		super().check_getindex_scalar(sigarray, refarray, index, result, refresult)
		assert isinstance(result, np.ndarray)
		assert result.dtype == sigarray.dtype

	def check_getindex_subseq(self, sigarray, refarray, index, result, refresult):
		super().check_getindex_subseq(sigarray, refarray, index, result, refresult)
		assert isinstance(result, SignatureArray)

	def check_getindex_slice(self, sigarray, refarray, index, result, refresult):
		super().check_getindex_slice(sigarray, refarray, index, result, refresult)



def test_uninitialized(refarray):
	"""Test creating with uninitialized() class method."""

	lengths = list(map(len, refarray))
	sigarray = SignatureArray.uninitialized(lengths)
	assert len(sigarray) == len(refarray)

	for i in range(len(sigarray)):
		assert sigarray.sizeof(i) == len(refarray[i])


def test_construct_from_signaturearray(sigarray):
	"""Test construction from another SignatureArray."""
	sa2 = SignatureArray(sigarray)
	assert sa2.dtype == sigarray.dtype
	assert sa2 == sigarray

	for dtype in map(np.dtype, ['i8', 'u4']):
		sa2 = SignatureArray(sigarray, dtype=dtype)
		assert sa2.dtype == dtype
		assert sa2 == sigarray


def test_empty():
	"""Really an edge case, but test it anyways."""

	sigarray = SignatureArray([])
	assert len(sigarray) == 0

	with pytest.raises(IndexError):
		sigarray[0]
