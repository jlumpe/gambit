"""Benchmarks for Jaccard distance metric."""

from pathlib import Path

import pytest
import numpy as np

from gambit.sigs import KmerSignature, load_signatures
from gambit.sigs.base import SignatureArray, AbstractSignatureArray
from gambit.sigs.hdf5 import HDF5Signatures
from gambit.metric import SCORE_DTYPE, jaccarddist_array
from gambit._cython.threads import omp_set_num_threads


def subset_sigs(sigs: AbstractSignatureArray, size: int, offset: int = 0) -> AbstractSignatureArray:
	length = len(sigs) - offset
	assert length >= size
	step = length // size
	subset = sigs[offset:offset + size * step:step]
	assert len(subset) == size
	return subset  # pyright: ignore[reportReturnType]


@pytest.fixture(scope='module', params=[100, 1000, 10000])
def refs_size(request):
	"""Number of reference signatures to process at a time.

	This corresponds to the ``chunksize`` parameter to :func:`gambit.metric.jaccarddist_matrix`,
	defaults to 1000 in :class:`gambit.query.QueryParams`.
	"""
	return request.param


@pytest.fixture(scope='module')
def ref_sigs_all(refseq_db_signatures_file: Path) -> HDF5Signatures:
	"""
	Full set of reference signatures RefSeq database (50,000 signatures), as HDF5Signatures object
	(not loaded into memory).
	"""
	return load_signatures(refseq_db_signatures_file)  # pyright: ignore[reportReturnType]


@pytest.fixture(scope='module')
def ref_sigs(ref_sigs_all: HDF5Signatures, refs_size: int) -> SignatureArray:
	"""Subset of reference signatures fully loaded into memory (contiguous Numpy array).

	Uses ``refs_size`` to determine the number of signatures to load into memory.
	"""
	return subset_sigs(ref_sigs_all, refs_size)  # pyright: ignore[reportReturnType]


@pytest.fixture(scope='module')
def query_sigs(ref_sigs_all: HDF5Signatures) -> list[KmerSignature]:
	"""Query signatures as individual Numpy arrays.

	Constant size, these will be used in a standard Python loop.
	"""
	# Use an offset so we're not taking the same ones as the references
	return list(subset_sigs(ref_sigs_all, 10, 10))


def _benchmark_jaccarddist_array(query_sigs: list[KmerSignature], ref_sigs: SignatureArray, out: np.ndarray):
	for query in query_sigs:
		jaccarddist_array(query, ref_sigs, out)


@pytest.mark.parametrize('threads', [1, 2, 4, 8])
def benchmark_jaccarddist_array(query_sigs: list[KmerSignature], ref_sigs: SignatureArray, benchmark, threads: int):
	"""Benchmark the jaccarddist_array function.

	This is the main workhorse function, ``jaccarddist_matrix()`` (used by ``query()``) calls this
	for each query signature against each chunk of reference signatures. The underlying Cython
	function should be parallelized.

	Use a preallocated array to store the results.
	"""

	omp_set_num_threads(threads)
	out = np.empty(len(ref_sigs), SCORE_DTYPE)
	benchmark(_benchmark_jaccarddist_array, query_sigs, ref_sigs, out)
