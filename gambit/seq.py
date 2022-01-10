"""Generic code for working with sequence data.

Note that all code in this package operates on DNA sequences as sequences of
bytes containing ascii-encoded nucleotide codes.

.. data:: NUCLEOTIDES

	``bytes`` corresponding to the four DNA nucleotides. Ascii-encoded upper
	case letters ``ACGT``. Note that the order, while arbitrary, is important
	in this variable as it defines how unique indices are assigned to k-mer
	sequences.
"""
from pathlib import Path
from typing import Union, Optional, IO, Iterable, List

from Bio import SeqIO
from Bio.Seq import Seq
from attr import attrs, attrib

from gambit._cython.kmers import revcomp
from gambit.util.io import FilePath
from gambit.util.io import open_compressed, ClosingIterator


# Byte representations of the four nucleotide codes in the order used for
# indexing k-mer sequences
NUCLEOTIDES = b'ACGT'

SEQ_TYPES = (str, bytes, bytearray, Seq)

#: Union of DNA sequence types accepted for k-mer search / signature calculation.
DNASeq = Union[SEQ_TYPES]

#: Sequence types accepted directly by native (Cython) code.
DNASeqBytes = Union[bytes, bytearray]


def seq_to_bytes(seq: DNASeq) -> DNASeqBytes:
	"""Convert generic DNA sequence to byte string representation.

	This is for passing sequence data to Cython functions.
	"""
	if isinstance(seq, (bytes, bytearray)):
		return seq
	if isinstance(seq, str):
		return seq.encode('ascii')
	if isinstance(seq, Seq):
		# This is recommended in the documentation over the deprecated encode() method, also
		# probably avoids copying any data as it typically just returns the seq._data attribute.
		return bytes(seq)
	raise TypeError(f'Expected sequence type, got {type(seq)}')


def validate_dna_seq_bytes(seq : bytes):
	"""Check that a sequence contains only valid nucleotide codes.

	Parameters
	----------
	seq : bytes
		ASCII-encoded nucleotide sequence.

	Raises
	------
	ValueError
		If the sequence contains an invalid nucleotide.
	"""
	for i, nuc in enumerate(seq):
		if nuc not in NUCLEOTIDES:
			raise ValueError(f'Invalid byte at position {i}: {nuc}')


@attrs(frozen=True, slots=True)
class SequenceFile:
	"""A reference to a DNA sequence file stored in the file system.

	Contains all the information needed to read and parse the file.

	Parameters
	----------
	path : Union[os.PathLike, str]
		Value of :attr:`path` attribute. May be string or path-like object.
	format : str
		Value of :attr:`format` attribute.
	compression : Optional[str]
		Value of :attr:`compression` attribute.

	Attributes
	----------
	path
		Path to the file.
	format
		String describing the file format as interpreted by
		:func:`Bio.SeqIO.parse`, e.g. ``'fasta'``.
	compression
		String describing compression method of the file, e.g. ``'gzip'``. None
		means no compression. See :func:`gambit.util.io.open_compressed`.
	"""
	path: Path = attrib(converter=Path)
	format: str = attrib()
	compression: Optional[str] = attrib(default=None)

	def open(self, mode: str = 'r', **kwargs) -> IO:
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
		IO
			Stream to file in given mode.
		"""
		return open_compressed(self.compression, self.path, mode, **kwargs)

	def parse(self, **kwargs) -> ClosingIterator[SeqIO.SeqRecord]:
		"""Open the file and lazily parse its contents.

		Returns iterator over sequence data in file. File is parsed lazily,
		and so must be kept open. The returned iterator is of type
		:class:`gambit.util.io.ClosingIterator` so it will close the file stream
		automatically when it finishes. It may also be used as a context manager
		that closes the stream on exit. You may also close the stream explicitly
		using the iterator's ``close`` method.

		Parameters
		----------
		\\**kwargs
			Keyword arguments to :meth:`open`.

		Returns
		-------
		gambit.util.io.ClosingIterator
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
	def from_paths(cls,
	               paths: Iterable[FilePath],
	               format: str,
	               compression: Optional[str] = None,
	               ) -> List['SequenceFile']:
		"""
		Create many instances at once from a collection of paths and a single
		format and compression type.

		Parameters
		----------
		paths
			Collection of paths as strings or path-like objects.
		format
			Sequence file format of files.
		compression
			Compression method of files.
		"""
		return [cls(path, format, compression) for path in paths]
