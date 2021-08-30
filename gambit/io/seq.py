"""Read and parse sequence files and calculate their k-mer signatures."""

from pathlib import Path
from typing import Optional, List, Sequence
from concurrent.futures import Executor, ThreadPoolExecutor, ProcessPoolExecutor, as_completed
from contextlib import nullcontext

import numpy as np
from Bio import SeqIO
from attr import attrs, attrib

from gambit.kmers import KmerSpec, find_kmers, dense_to_sparse, KmerSignature
from .util import open_compressed, ClosingIterator
from gambit.util.progress import iter_progress, get_progress


@attrs(frozen=True, slots=True)
class SequenceFile:
	"""A reference to a DNA sequence file stored in the file system.

	Contains all the information needed to read and parse the file.

	Parameters
	----------
	path : Union[pathlib.Path, str]
		Value of :attr:`path` attribute. May be string or path-like object.
	format : str
		Value of :attr:`format` attribute.
	compression : Optional[str]
		Value of :attr:`compression` attribute.

	Attributes
	----------
	path : pathlib.Path
		Path to the file.
	format : str
		String describing the file format as interpreted by
		:func:`Bio.SeqIO.parse`, e.g. ``'fasta'``.
	compression : str or None
		String describing compression method of the file, e.g. ``'gzip'``. None
		means no compression. See :func:`gambit.io.util.open_compressed`.
	"""
	path: Path = attrib(converter=Path)
	format: str = attrib()
	compression: Optional[str] = attrib(default=None)

	def open(self, mode: str = 'r', **kwargs):
		"""
		Open a stream to the file, with compression/decompression applied
		transparently.

		Parameters
		----------

		mode : str
			Same as equivalent argument to the built-in :func:open`. Some modes may not be supported
			by all compression types.
		\\**kwargs
			Additional text mode specific keyword arguments to pass to opener. Equivalent to the
			following arguments of the built-in :func:`open`: ``encoding``, ``errors``, and
			``newlines``. May not be supported by all compression types.

		Returns
		-------
		Stream to file in given mode.
		"""
		return open_compressed(self.compression, self.path, mode, **kwargs)

	def parse(self, **kwargs):
		"""Open the file and lazily parse its contents.

		Returns iterator over sequence data in file. File is parsed lazily,
		and so must be kept open. The returned iterator is of type
		:class:`gambit.io.util.ClosingIterator` so it will close the file stream
		automatically when it finishes. It may also be used as a context manager
		that closes the stream on exit. You may also close the stream explicitly
		using the iterator's ``close`` method.

		Parameters
		----------
		\\**kwargs
			Keyword arguments to :meth:`open`.

		Returns
		-------
		gambit.io.util.ClosingIterator
			Iterator yielding :class:`Bio.SeqIO.SeqRecord` instances for each sequence in the file.
		"""

		fobj = self.open('rt', **kwargs)

		try:
			records = SeqIO.parse(fobj, self.format)
			return ClosingIterator(records, fobj)

		except:
			fobj.close()
			raise

	def absolute(self) -> 'SequenceFile':
		"""Make a copy of the instance with an absolute path."""
		if self.path.is_absolute():
			return self
		else:
			return SequenceFile(self.path.absolute(), self.format, self.compression)

	@classmethod
	def from_paths(cls, paths, format: str, compression: Optional[str] = None) -> List['SequenceFile']:
		"""
		Create many instances at once from a collection of paths and a single
		format and compression type.

		Parameters
		----------
		paths
			Collection of paths as strings or path-like objects.
		format : str
			Sequence file format of files.
		compression : str
			Compression method of files.
		"""
		return [cls(path, format, compression) for path in paths]


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
		Pre-allocated numpy array to write dense output to. Should be of length ``kspec.idx_len``.
		Note that this is still used as working space even if ``sparse=True``. Should be zeroed
		prior to use (although if not the result will effectively be the bitwise AND between its
		previous value and k-mers found in ``data``.

	Returns
	-------
	numpy.ndarray
		If ``sparse`` is False, returns dense K-mer vector (same array as ``dense_out`` if it was
		given). If ``sparse`` is True returns k-mers in sparse coordinate format (dtype will match
		:func:`gambit.kmers.vec_to_coords`).

	See Also
	--------
	gambit.kmers.find_kmers
	.find_kmers_in_file
	"""
	if dense_out is None:
		dense_out = np.zeros(kspec.idx_len, dtype=bool)

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
		:func:`gambit.kmers.vec_to_coords`).

	See Also
	--------
	gambit.kmers.find_kmers
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
