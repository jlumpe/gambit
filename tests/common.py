"""Helper functions for tests."""

from typing import Optional, Union, Sequence

import numpy as np

from gambit.kmers import KmerSpec, kmer_to_index
from gambit.seq import seq_to_bytes, revcomp
from gambit.sigs import KmerSignature, SignatureArray
from gambit.sigs.calc import dense_to_sparse, sparse_to_dense
from gambit.db import Taxon


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


def make_signatures(k_or_kspec: Union[int, KmerSpec], n: int, dtype: np.dtype = np.dtype('u8')) -> SignatureArray:
	"""Make artificial k-mer signatures.

	Parameters
	----------
	k_or_kspec
	n
		Number of signatures to create.
	dtype
		Numpy dtype of signatures.
	"""
	if isinstance(k_or_kspec, int):
		kspec = KmerSpec(k_or_kspec, 'A')
	elif isinstance(k_or_kspec, KmerSpec):
		kspec = k_or_kspec
	else:
		raise TypeError(f'Expected int or KmerSpec, got {type(k_or_kspec)}')

	nk = kspec.nkmers

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

	return SignatureArray(signatures_list, kspec, dtype=dtype)


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


def make_kmer_seq(kspec: KmerSpec,
                  seqlen: int,
                  kmer_interval: int,
                  n_interval: Optional[int] = None,
                  ) -> tuple[bytes, KmerSignature]:
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


def make_kmer_seqs(kspec: KmerSpec,
                   nseqs: int,
                   seqlen: int,
                   kmer_interval: int,
                   n_interval: Optional[int] = None,
                   ) -> tuple[list[bytes], KmerSignature]:
	"""Create a set of DNA sequences with known combined signature."""

	seqs = []
	vec = np.zeros(kspec.nkmers, dtype=bool)

	for i in range(nseqs):
		seq, sig = make_kmer_seq(kspec, seqlen, kmer_interval, n_interval)

		# Combine vectors of all sequences
		vec |= sparse_to_dense(kspec, sig)

		# Convert every other sequence to lower case, just to switch things up...
		if i % 2:
			seq = seq.lower()

		seqs.append(seq)

	return seqs, dense_to_sparse(vec)


def make_lineage(thresholds: Sequence[Optional[float]]) -> list[Taxon]:
	"""Create a linage of taxa that have the given distance thresholds.

	Parameters
	----------
	thresholds
		Distance thresholds for taxa, in ascending order.

	Returns
	-------
	List[Taxon]
		Created taxa in ascending order.
	"""
	taxa = []

	for i, t in enumerate(thresholds):
		taxon = Taxon(
			name=f't{i + 1}',
			distance_threshold=t,
		)
		if taxa:
			taxa[-1].parent = taxon
		taxa.append(taxon)

	return taxa
