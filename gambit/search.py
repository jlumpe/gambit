"""Search for k-mer matches in sequence data."""

from typing import Union, Optional, List, Sequence
from concurrent.futures import Executor, ThreadPoolExecutor, ProcessPoolExecutor, as_completed
from contextlib import nullcontext

import numpy as np
from Bio import SeqIO

from gambit.kmers import NUCLEOTIDES, KmerSpec, KmerSignature, dense_to_sparse, kmer_to_index, reverse_complement
from gambit.io.seq import SequenceFile
from gambit.util.progress import iter_progress, get_progress


def find_kmers(kspec: KmerSpec,
               seq: Union[bytes, str],
               *,
               sparse: bool = True,
               dense_out: Optional[np.ndarray] = None,
               ) -> KmerSignature:
	"""Find all k-mers in a DNA sequence.

	Searches sequence both backwards and forwards (reverse complement). The sequence may contain
	invalid characters (not one of the four nucleotide codes) which will simply not be matched.

	Parameters
	----------
	kspec : .KmerSpec
		K-mer spec to use for search.
	seq
		Sequence to search within as ``bytes`` or ``str``. If ``str`` will be encoded as ASCII.
		Lower-case characters are OK and will be matched as upper-case.
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
	.find_kmers_parse
	"""
	if dense_out is None:
		dense_out = np.zeros(kspec.nkmers, dtype=bool)

	# Convert sequence to bytes
	if not isinstance(seq, bytes):
		if not isinstance(seq, str):
			seq = str(seq)

		seq = seq.encode('ascii')

	# Convert to upper-case only if needed
	nucs_lower = set(NUCLEOTIDES.lower())
	for char in seq:
		if char in nucs_lower:
			seq = seq.upper()
			break

	_find_kmers(kspec, seq, dense_out)

	if sparse:
		return dense_to_sparse(dense_out)
	else:
		return dense_out


def _find_kmers(kspec, seq, out):
	"""Actual implementation of find_kmers.

	Parameters
	----------
	kspec : KmerSpec
	seq : bytes
		Upper-case ASCII nucleotide codes.
	out : np.ndarray
		Write dense output to this array.
	"""

	# Reverse complement of prefix
	rcprefix = reverse_complement(kspec.prefix)

	# Search forward
	start = 0
	while True:
		loc = seq.find(kspec.prefix, start, -kspec.k)
		if loc < 0:
			break

		kmer = seq[loc + kspec.prefix_len:loc + kspec.total_len]
		if not isinstance(kmer, bytes):
			kmer = str(kmer).encode('ascii')

		try:
			out[kmer_to_index(kmer)] = 1
		except ValueError:
			pass

		start = loc + 1

	# Search backward
	start = kspec.k
	while True:
		loc = seq.find(rcprefix, start)
		if loc < 0:
			break

		rckmer = seq[loc - kspec.k:loc]
		if not isinstance(rckmer, bytes):
			rckmer = str(rckmer).encode('ascii')
		kmer = reverse_complement(rckmer)

		try:
			out[kmer_to_index(kmer)] = 1
		except ValueError:
			pass

		start = loc + 1


def find_kmers_parse(kspec: KmerSpec, data, format: str, *, sparse: bool = True, dense_out: Optional[np.ndarray] = None) -> np.ndarray:
	"""Parse sequence data with ``Bio.Seq.parse()`` and find k-mers.

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
	.find_kmers
	.find_kmers_in_file
	"""
	if dense_out is None:
		dense_out = np.zeros(kspec.nkmers, dtype=bool)

	for record in SeqIO.parse(data, format):
		find_kmers(kspec, record.seq, dense_out=dense_out)

	if sparse:
		return dense_to_sparse(dense_out)
	else:
		return dense_out


def find_kmers_in_file(kspec: KmerSpec, seqfile: SequenceFile, *, sparse: bool = True, dense_out: Optional[np.ndarray] = None) -> np.ndarray:
	"""Open a sequence file on disk and find k-mers.

	This works identically to :func:`.find_kmers_parse` but takes a :class:`.SequenceFile` as input
	instead of a data stream.

	Parameters
	----------
	kspec
		Spec for k-mer search.
	seqfile
		File to read.
	sparse
		See :func:`.find_kmers_parse`.
	dense_out
		See :func:`.find_kmers_parse`.

	Returns
	-------
	numpy.ndarray
		If ``sparse`` is False, returns dense K-mer vector (same array as ``dense_out`` if it was
		given). If ``sparse`` is True returns k-mers in sparse coordinate format (dtype will match
		:func:`gambit.kmers.dense_to_sparse`).

	See Also
	--------
	.find_kmers
	.find_kmers_in_files
	.find_kmers_parse
	"""
	with seqfile.open() as f:
		return find_kmers_parse(kspec, f, seqfile.format, sparse=sparse, dense_out=dense_out)


def find_kmers_in_files(kspec: KmerSpec,
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
	.find_kmers_in_file
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
				sigs.append(find_kmers_in_file(kspec, file))

	else:
		sigs = [None] * len(files)
		future_to_index = dict()

		with executor_context, get_progress(progress, len(files)) as meter:
			for i, file in enumerate(files):
				future = executor.submit(find_kmers_in_file, kspec, file)
				future_to_index[future] = i

			for future in as_completed(future_to_index):
				i = future_to_index[future]
				sigs[i] = future.result()
				meter.increment()

		assert all(sig is not None for sig in sigs)

	return sigs
