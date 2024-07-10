"""Calculate k-mer signatures from sequence data."""

from typing import Optional, Sequence, MutableSet, Union, Iterable
from abc import abstractmethod
from concurrent.futures import Executor, ThreadPoolExecutor, ProcessPoolExecutor, as_completed
from contextlib import nullcontext

import numpy as np

from .base import KmerSignature, SignatureList
from gambit.kmers import KmerSpec, find_kmers, kmer_to_index, nkmers, index_dtype
from gambit.seq import SEQ_TYPES, DNASeq, SequenceFile
from gambit.util.progress import iter_progress, get_progress


class KmerAccumulator(MutableSet[int]):
	"""Base class for data structures which track k-mers as they are found in sequences.

	Implements the ``MutableSet`` interface for k-mer indices. Indices are added via :meth:`add` or
	:meth:`add_kmer` methods, when finished a sparse k-mer signature can be obtained from
	:meth:`signature`.
	"""
	k: int

	def add_kmer(self, kmer: bytes):
		"""Add a k-mer by its sequence rather than its index.

		Argument may contain invalid (non-nucleotide) bytes, in which case it is ignored.
		"""
		if len(kmer) != self.k:
			raise ValueError(f'Expected {self.k}-mer, argument has length {len(kmer)}')

		try:
			idx = kmer_to_index(kmer)
		except ValueError:
			return

		self.add(idx)

	@abstractmethod
	def signature(self) -> KmerSignature:
		"""Get signature for accumulated k-mers."""
		pass


class ArrayAccumulator(KmerAccumulator):
	"""K-mer accumulator implemented as a dense boolean array.

	This is pretty efficient for smaller values of ``k``, but time and space requirements increase
	exponentially with larger values.
	"""
	array: np.ndarray

	def __init__(self, k: int):
		self.k = k
		self.array = np.zeros(nkmers(k), dtype=bool)
		self._dtype = index_dtype(self.k)

	def __len__(self):
		return self.array.sum()

	def __iter__(self):
		for i, x in enumerate(self.array):
			if x:
				yield self._dtype.type(i)

	def __contains__(self, index: int):
		return self.array[index]

	def add(self, i: int):
		self.array[i] = True

	def discard(self, i: int):
		self.array[i] = False

	def clear(self):
		self.array[:] = False

	def signature(self) -> KmerSignature:
		return np.flatnonzero(self.array).astype(self._dtype)


class SetAccumulator(KmerAccumulator):
	"""Accumulator which uses the builtin Python ``set`` class.

	This has more overhead than the array version for smaller values of ``k`` but behaves much
	better asymptotically.
	"""

	set: set

	def __init__(self, k: int):
		self.k = k
		self.set = set()
		self._dtype = index_dtype(self.k)

	def __len__(self):
		return len(self.set)

	def __iter__(self):
		return iter(self.set)

	def __contains__(self, index):
		return index in self.set

	def clear(self):
		self.set.clear()

	def discard(self, index: int):
		self.set.discard(index)

	def add(self, index: int):
		self.set.add(self._dtype.type(index))

	def signature(self) -> KmerSignature:
		sig = np.fromiter(self.set, dtype=self._dtype)
		sig.sort()
		return sig


def default_accumulator(k: int) -> KmerAccumulator:
	"""Get a default k-mer accumulator instance for the given value of ``k``.

	Returns a :class:`.ArrayAccumulator` for ``k <= 11`` and a :class:`.SetAccumulator` for
	``k > 11``.
	"""
	return SetAccumulator(k) if k > 11 else ArrayAccumulator(k)


def accumulate_kmers(accumulator: KmerAccumulator, kmerspec: KmerSpec, seq: DNASeq):
	"""Find k-mer matches in sequence and add their indices to an accumulator."""
	for match in find_kmers(kmerspec, seq):
		try:
			index = match.kmer_index()
		except ValueError:
			continue
		accumulator.add(index)


def calc_signature(kmerspec: KmerSpec,
                   seqs: Union[DNASeq, Iterable[DNASeq]],
                   *,
                   accumulator: Optional[KmerAccumulator] = None,
                   ) -> KmerSignature:
	"""Calculate the k-mer signature of a DNA sequence or set of sequences.

	Searches sequences both backwards and forwards (reverse complement). Sequences may contain
	invalid characters (not one of the four nucleotide codes) which will simply not be matched.

	Parameters
	----------
	kmerspec
		K-mer spec to use for search.
	seqs
		Sequence or sequences to search within. Lowercase characters are OK.
	accumulator
		TODO

	Returns
	-------
	numpy.ndarray
		K-mer signature in sparse coordinate format. Data type will be ``kspec.index_dtype``.

	See Also
	--------
	.calc_file_signature
	"""
	if isinstance(seqs, SEQ_TYPES):
		seqs = [seqs]

	if accumulator is None:
		accumulator = default_accumulator(kmerspec.k)

	for seq in seqs:
		accumulate_kmers(accumulator, kmerspec, seq)

	return accumulator.signature()


def calc_file_signature(kspec: KmerSpec,
                        seqfile: SequenceFile,
                        *,
                        accumulator: Optional[KmerAccumulator] = None,
                        ) -> KmerSignature:
	"""Open a sequence file on disk and calculate its k-mer signature.

	This works identically to :func:`.calc_signature_parse` but takes a :class:`.SequenceFile` as
	input instead of a data stream.

	Parameters
	----------
	kspec
		Spec for k-mer search.
	seqfile
		File to read.
	accumulator
		TODO

	Returns
	-------
	numpy.ndarray
		K-mer signature in sparse coordinate format (dtype will match
		:func:`gambit.kmers.dense_to_sparse`).

	See Also
	--------
	.calc_signature
	.calc_file_signatures
	"""
	with seqfile.parse() as records:
		return calc_signature(kspec, (record.seq for record in records))


def calc_file_signatures(kspec: KmerSpec,
                         files: Sequence[SequenceFile],
                         progress=None,
                         concurrency: Optional[str] = 'processes',
                         max_workers: Optional[int] = None,
                         executor: Optional[Executor] = None,
                         ) -> SignatureList:
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
		Process files concurrently. ``"processes"`` for process-based (default), ``"threads"`` for
		threads-based, ``None`` for no concurrency.
	max_workers
		Number of worker threads/processes to use if ``concurrency`` is not None.
	executor
		Instance of class:`concurrent.futures.Executor` to use for concurrency. Overrides the
		``concurrency`` and ``max_workers`` arguments.

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

	return SignatureList(sigs, kspec)
