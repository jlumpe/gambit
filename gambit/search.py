"""Search for k-mer matches in sequence data."""

from typing import Optional, List, Sequence, Iterator
from concurrent.futures import Executor, ThreadPoolExecutor, ProcessPoolExecutor, as_completed
from contextlib import nullcontext

import numpy as np
from Bio import SeqIO
from attr import attrs, attrib

from gambit.kmers import NUCLEOTIDES, DNASeq, KmerSpec, KmerSignature, dense_to_sparse, \
	kmer_to_index, kmer_to_index_rc, seq_to_bytes, revcomp
from gambit.io.seq import SequenceFile
from gambit.util.progress import iter_progress, get_progress


@attrs(slots=True)
class KmerMatch:
	"""Represents a

	Attributes
	----------
	kmerspec
		K-mer spec used for search.
	seq
		The sequence searched within.
	pos
		Index of first nucleotide of prefix in ``seq``.
	reverse
		If the match is on the reverse strand.
	"""
	kmerspec: KmerSpec = attrib()
	seq: DNASeq = attrib()
	pos: int = attrib()
	reverse: bool = attrib()

	def kmer_indices(self) -> slice:
		"""Index range for k-mer in sequence (without prefix)."""
		if self.reverse:
			return slice(self.pos - self.kmerspec.total_len + 1, self.pos - self.kmerspec.prefix_len + 1)
		else:
			return slice(self.pos + self.kmerspec.prefix_len, self.pos + self.kmerspec.total_len)

	def full_indices(self) -> slice:
		"""Index range for prefix plus k-mer in sequence."""
		if self.reverse:
			return slice(self.pos - self.kmerspec.total_len + 1, self.pos + 1)
		else:
			return slice(self.pos, self.pos + self.kmerspec.total_len)

	def kmer(self) -> bytes:
		"""Get matched k-mer sequence."""
		kmer = seq_to_bytes(self.seq[self.kmer_indices()])
		return revcomp(kmer) if self.reverse else kmer

	def kmer_index(self) -> int:
		"""Get index of matched k-mer.

		Raises
		------
		ValueError
			If the k-mer contains invalid nucleotides.
		"""
		kmer = self.seq[self.kmer_indices()]
		return kmer_to_index_rc(kmer) if self.reverse else kmer_to_index(kmer)


def find_kmers(kmerspec: KmerSpec, seq: DNASeq) -> Iterator[KmerMatch]:
	"""Locate k-mers with the given prefix in a DNA sequence.

	Searches sequence both backwards and forwards (reverse complement). The sequence may contain
	invalid characters (not one of the four nucleotide codes) which will simply not be matched.

	Parameters
	----------
	kmerspec
		K-mer spec to use for search.
	seq
		Sequence to search within. Lowercase characters are OK and will be matched as uppercase.

	Returns
	-------
	Iterator[KmerMatch]
		Iterator of :class:`.KmerMatch` objects.
	"""

	haystack = seq_to_bytes(seq)

	# Convert to uppercase only if needed
	nucs_lower = NUCLEOTIDES.lower()
	for char in haystack:
		if char in nucs_lower:
			haystack = haystack.upper()
			break

	# Find forward
	start = 0

	while True:
		loc = haystack.find(kmerspec.prefix, start, -kmerspec.k)
		if loc < 0:
			break

		yield KmerMatch(kmerspec, seq, loc, False)

		start = loc + 1

	# Find reverse
	prefix_rc = revcomp(kmerspec.prefix)
	start = kmerspec.k

	while True:
		loc = haystack.find(prefix_rc, start)
		if loc < 0:
			break

		yield KmerMatch(kmerspec, seq, loc + kmerspec.prefix_len - 1, True)

		start = loc + 1


def calc_signature(kmerspec: KmerSpec,
                   seq: DNASeq,
                   *,
                   sparse: bool = True,
                   dense_out: Optional[np.ndarray] = None,
                   ) -> KmerSignature:
	"""Calculate the k-mer signature of a DNA sequence.

	Searches sequence both backwards and forwards (reverse complement). The sequence may contain
	invalid characters (not one of the four nucleotide codes) which will simply not be matched.

	Parameters
	----------
	kmerspec : .KmerSpec
		K-mer spec to use for search.
	seq
		Sequence to search within. Lowercase characters are OK and will be matched as uppercase.
	dense_out : numpy.ndarray
		Pre-allocated numpy array to write dense output to. Should be of length ``kspec.nkmers``.
		Note that this is still used as working space even if ``sparse=True``. Should be zeroed
		prior to use (although if not the result will effectively be the bitwise AND between its
		previous value and k-mers found in ``data``.
	sparse : bool
		If True return k-mers in sparse coordinate format rather than dense (bit vector) format.

	Returns
	-------
	numpy.ndarray
		If ``sparse`` is False, returns dense K-mer vector (same array as ``dense_out`` if it was
		given). If ``sparse`` is True returns k-mers in sparse coordinate format (dtype will match
		:func:`gambit.kmers.dense_to_sparse`).

	See Also
	--------
	.calc_signature_parse
	"""
	if dense_out is None:
		dense_out = np.zeros(kmerspec.nkmers, dtype=bool)

	for match in find_kmers(kmerspec, seq):
		try:
			i = match.kmer_index()
		except ValueError:
			continue
		dense_out[i] = 1

	if sparse:
		return dense_to_sparse(dense_out)
	else:
		return dense_out


def calc_signature_parse(kspec: KmerSpec, data, format: str, *, sparse: bool = True, dense_out: Optional[np.ndarray] = None) -> np.ndarray:
	"""Parse sequence data with ``Bio.Seq.parse()`` and calculate its k-mer signature.

	Parameters
	----------
	kspec : gambit.kmers.KmerSpec
		Spec for k-mer search.
	data
		Stream with sequence data. Readable file-like object in text mode.
	format : str
		Sequence file format, as interpreted by :func:`Bio.SeqIO.parse`.
	sparse : bool
		If True return k-mers in sparse coordinate format rather than dense (bit vector) format.
	dense_out : numpy.ndarray
		Pre-allocated numpy array to write dense output to. Should be of length ``kspec.nkmers``.
		Note that this is still used as working space even if ``sparse=True``. Should be zeroed
		prior to use (although if not the result will effectively be the bitwise AND between its
		previous value and k-mers found in ``data``.

	Returns
	-------
	numpy.ndarray
		If ``sparse`` is False, returns dense K-mer vector (same array as ``dense_out`` if it was
		given). If ``sparse`` is True returns k-mers in sparse coordinate format (dtype will match
		:func:`gambit.kmers.dense_to_sparse`).

	See Also
	--------
	.calc_signature
	.calc_file_signature
	"""
	if dense_out is None:
		dense_out = np.zeros(kspec.nkmers, dtype=bool)

	for record in SeqIO.parse(data, format):
		calc_signature(kspec, record.seq, dense_out=dense_out)

	if sparse:
		return dense_to_sparse(dense_out)
	else:
		return dense_out


def calc_file_signature(kspec: KmerSpec, seqfile: SequenceFile, *, sparse: bool = True, dense_out: Optional[np.ndarray] = None) -> np.ndarray:
	"""Open a sequence file on disk and calculate its k-mer signature.

	This works identically to :func:`.calc_signature_parse` but takes a :class:`.SequenceFile` as
	input instead of a data stream.

	Parameters
	----------
	kspec
		Spec for k-mer search.
	seqfile
		File to read.
	sparse
		See :func:`.calc_signature_parse`.
	dense_out
		See :func:`.calc_signature_parse`.

	Returns
	-------
	numpy.ndarray
		If ``sparse`` is False, returns dense K-mer vector (same array as ``dense_out`` if it was
		given). If ``sparse`` is True returns k-mers in sparse coordinate format (dtype will match
		:func:`gambit.kmers.dense_to_sparse`).

	See Also
	--------
	.calc_signature
	.calc_file_signatures
	.calc_signature_parse
	"""
	with seqfile.open() as f:
		return calc_signature_parse(kspec, f, seqfile.format, sparse=sparse, dense_out=dense_out)


def calc_file_signatures(kspec: KmerSpec,
                         files: Sequence[SequenceFile],
                         progress=None,
                         concurrency: Optional[str] = 'threads',
                         max_workers: Optional[int] = None,
                         executor: Optional[Executor] = None,
                         ) -> List[KmerSignature]:
	"""Parse and calculate k-mer signatures for multiple sequence files.

	Parameters
	----------
	kspec
		Spec for k-mer search.
	seqfile
		Files to read.
	progress
		Display a progress meter. See :func:`gambit.util.progress.get_progress` for allowed values.
	concurrency
		Process files concurrently. ``"threads"`` for thread-based (default), ``"processes"`` for
		process-based, ``None`` for no concurrency.
	max_workers
		Number of worker threads/processes to use if ``concurrency`` is not None.
	executor
		Instance of class:`concurrent.futures.Executor` to use for concurrency. Overrides the
		``concurrency`` and ``max_workers`` arguments.

	Returns
	-------
		List of signatures in sparse coordinate format.

	See Also
	--------
	.calc_file_signature
	"""
	if executor is None:
		if concurrency == 'threads':
			executor = ThreadPoolExecutor(max_workers=max_workers)
		elif concurrency == 'processes':
			executor = ProcessPoolExecutor(max_workers=max_workers)
		elif concurrency is not None:
			raise ValueError(f'concurrency should be one of [None, "threads", "processes"], got {concurrency!r}')

		executor_context = executor

	else:
		executor_context = nullcontext()

	if executor is None:
		sigs = []

		with iter_progress(files, progress) as file_itr:
			for file in file_itr:
				sigs.append(calc_file_signature(kspec, file))

	else:
		sigs = [None] * len(files)
		future_to_index = dict()

		with executor_context, get_progress(progress, len(files)) as meter:
			for i, file in enumerate(files):
				future = executor.submit(calc_file_signature, kspec, file)
				future_to_index[future] = i

			for future in as_completed(future_to_index):
				i = future_to_index[future]
				sigs[i] = future.result()
				meter.increment()

		assert all(sig is not None for sig in sigs)

	return sigs
