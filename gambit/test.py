"""Helper functions for tests."""

from typing import Optional, Tuple, Union, ContextManager
from contextlib import contextmanager

import numpy as np
from Bio.Seq import Seq

from gambit.kmers import KmerSpec, KmerSignature, dense_to_sparse, kmer_to_index, revcomp, nkmers, \
	seq_to_bytes
from gambit.signatures import SignatureArray
from gambit.query import QueryResultItem
from gambit.classify import ClassifierResult, GenomeMatch
from gambit.util.progress import TestProgressMeter, ProgressConfig, progress_config, capture_progress


# Sequence types used for k-mer search
SEQ_TYPES = [str, bytes, bytearray, Seq]

def convert_seq(seq, type):
	"""Convert sequence to any of the accepted argument types for k-mer search."""
	seq = seq_to_bytes(seq)
	if type is str:
		return seq.decode('ascii')
	return type(seq)


def bernoulli(size: Union[int, tuple], p: float) -> np.ndarray:
	"""Sample from Bernoulli distribution using Numpy.

	Parameters
	----------
	size
		Size of output array.
	p
		Probability of True.
	"""
	return np.random.choice([False, True], size, p=[1 - p, p])


def make_signatures(k: int, n: int, dtype: np.dtype = np.dtype('u8')) -> SignatureArray:
	"""Make artificial k-mer signatures.

	Parameters
	----------
	k
	n
		Number of signatures to create.
	dtype
		Numpy dtype of signatures.
	"""
	nk = nkmers(k)

	# Uniform probability for including a k-mer
	# .005 is about what we see in real genomes with k=11
	# For small values of k, try to have an expected signature length of at least 20 but don't go over p=.5
	p = max(.005, min(20 / nk, .5))

	signatures_list = []

	# Add empty and full sets as edge cases
	signatures_list.append(np.arange(0))
	signatures_list.append(np.arange(nk))

	# Use a core set of k-mers so that we get some overlap
	core_vec = bernoulli(nk, p)

	for i in range(n - 3):
		vec = bernoulli(nk, p) | core_vec
		signatures_list.append(dense_to_sparse(vec))

	# Add one more that does not include core set
	vec = bernoulli(nk, p) & ~core_vec
	signatures_list.append(dense_to_sparse(vec))

	return SignatureArray(signatures_list, dtype=dtype)


def random_seq(n: int, chars: str = 'ACGT') -> bytes:
	"""Generate a simple random DNA sequence.

	Parameters
	----------
	n : int
		Length of sequence to generate.
	chars : str
		Characters to use for sequence. Must be encodable as ascii.

	Returns
	-------
		Sequence as bytes or str.
	"""
	chars_array = np.frombuffer(chars.encode('ascii'), dtype='u1')
	return np.random.choice(chars_array, n).tobytes()


def fill_bytearray(pattern: bytes, n: int) -> bytearray:
	"""Create a bytearray with a repeating pattern.

	Parameters
	----------
	pattern : bytes
		Pattern to repeat in array.
	n : int
		Length of array to create (not number of repeats).

	Returns
	-------
	bytearray
		Filled array
	"""
	array = bytearray(n)
	p = len(pattern)

	for i in range(n):
		array[i] = pattern[i % p]

	return array


def make_kmer_seq(kspec: KmerSpec, seqlen: int, kmer_interval: int, n_interval: Optional[int] = None
                  ) -> Tuple[bytes, KmerSignature]:
	"""Create a DNA sequence with a known k-mer signature.

	The sequence consists of a background of N's with a k-mer match every ``kmer_interval``
	nucleotides.

	Parameters
	----------
	kspec
	seqlen
		Length of sequence.
	kmer_interval
		Number of nucleotides between each k-mer match added.
	n_interval
		Every this many k-mers, add an N to the k-mer sequence to create a k-mer that should not be
		matched.

	Returns
	-------
	tuple
		Tuple of (seq, signature).
	"""
	if kmer_interval < kspec.total_len:
		raise ValueError()

	# Initialize filled with N's
	seq_array = fill_bytearray(b'N', seqlen)

	# Keep track of which kmers have been added
	vec = np.zeros(kspec.nkmers, dtype=bool)

	# Add matches
	for i, p in enumerate(range(0, seqlen - kspec.total_len, kmer_interval)):
		# Pick random k-mer, but make sure its reverse complement doesn't cause another match.
		while True:
			kmer = random_seq(kspec.k)
			if not kmer.endswith(revcomp(kspec.prefix)):
				break

		# Every so often add an N just to throw things off
		if n_interval is not None and i % n_interval == 0:
			kmer_array = bytearray(kmer)
			kmer_array[np.random.randint(kspec.k)] = ord(b'N')
			kmer = bytes(kmer_array)

		else:
			vec[kmer_to_index(kmer)] = True

		match = kspec.prefix + kmer

		# Reverse every other match
		if i % 2 == 1:
			match = revcomp(match)

		seq_array[p:p + kspec.total_len] = match

	return bytes(seq_array), dense_to_sparse(vec)


def compare_genome_matches(match1: Optional[GenomeMatch], match2: Optional[GenomeMatch]) -> bool:
	"""Compare two ``GenomeMatch`` instances for equality.

	The values for the ``distance`` attribute are only checked for approximate equality, to support
	instances where one was loaded from a results archive (saving and loading a float in JSON is
	lossy).

	Also allows one or both values to be None.
	"""
	if match1 is None or match2 is None:
		return match1 is None and match2 is None

	return match1.genome == match2.genome and \
	       match1.matched_taxon == match2.matched_taxon and \
	       np.isclose(match1.distance, match2.distance)


def compare_classifier_results(result1: ClassifierResult, result2: ClassifierResult) -> bool:
	"""Compare two ``ClassifierResult`` instances for equality."""
	return result1.success == result2.success and \
	       result1.predicted_taxon == result2.predicted_taxon and \
	       compare_genome_matches(result1.primary_match, result2.primary_match) and \
	       compare_genome_matches(result1.closest_match, result2.closest_match) and \
	       set(result1.warnings) == set(result2.warnings) and \
	       result1.error == result2.error


def compare_result_items(item1: QueryResultItem, item2: QueryResultItem) -> bool:
	"""Compare two ``QueryResultItem`` instances for equality.

	Does not compare the value of the ``input`` attributes.
	"""
	return item1.report_taxon == item2.report_taxon and \
	       compare_classifier_results(item1.classifier_result, item2.classifier_result)


@contextmanager
def check_progress(*,
                   total: Optional[int] = None,
                   allow_decrement: bool = False,
                   check_closed: bool = True,
                   ) -> ContextManager[ProgressConfig]:
	"""Context manager which checks a progress meter is advanced to completion.

	Returned context manager yields a ``ProgressConfig`` instance on enter, tests are run when
	context is exited. Expects that the config will be used to instantiate exactly one progress
	meter. Tests are performed with assert statements.

	Parameters
	----------
	total
		Check that the progress meter is created with this total length.
	allow_decrement
		If false, raise an error if the created progress meter is moved backwards.
	check_closed
		Check that the progress meter was closed.
	"""

	conf = progress_config(TestProgressMeter, allow_decrement=allow_decrement)
	conf2, l = capture_progress(conf)

	yield conf2

	assert len(l) != 0, 'Progress meter not instantiated'
	assert len(l) == 1, 'Progress meter instantiated multiple times'

	pbar = l[0]

	assert pbar.n == pbar.total, 'Progress meter not completed'

	if total is not None:
		assert pbar.total == total

	if check_closed:
		assert pbar.closed
