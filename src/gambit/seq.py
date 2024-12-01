"""Generic code for working with sequence data.

Note that all code in this package operates on DNA sequences as sequences of
bytes containing ascii-encoded nucleotide codes.


.. data:: NUCLEOTIDES

	``bytes`` corresponding to the four DNA nucleotides. Ascii-encoded upper
	case letters ``ACGT``. Note that the order, while arbitrary, is important
	in this variable as it defines how unique indices are assigned to k-mer
	sequences.

.. class:: DNASeq

	Type alias for DNA sequence types accepted for k-mer search / signature calculation
	(``str``, ``bytes``, ``bytearray``, or :class:`Bio.Seq.Seq`).
"""

from pathlib import Path
from typing import Union, Optional, IO, Iterable
from os import PathLike

from Bio import SeqIO
from Bio.Seq import Seq
from attr import attrs, attrib
from typing_extensions import TypeAlias

from gambit._cython.kmers import revcomp
from gambit.util.io import FilePath
from gambit.util.io import open_compressed, ClosingIterator


# Byte representations of the four nucleotide codes in the order used for
# indexing k-mer sequences
NUCLEOTIDES = b'ACGT'

SEQ_TYPES = (str, bytes, bytearray, Seq)

DNASeq: TypeAlias = Union[SEQ_TYPES]
# Type alias for sequence types accepted directly by native (Cython) code.
DNASeqBytes: TypeAlias = Union[bytes, bytearray]


def seq_to_bytes(seq: 'DNASeq') -> 'DNASeqBytes':
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


def validate_dna_seq_bytes(seq: DNASeqBytes):
	"""Check that a sequence contains only valid nucleotide codes (upper case).

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


def parse_seqs(path: FilePath,
               format: str = 'fasta',
               compression: str = 'auto',
               **kwargs) -> ClosingIterator[SeqIO.SeqRecord]:
	"""Open a sequence file and lazily parse its contents.

	This is essentially a wrapper over BioPython's :func:`Bio.SeqIO.parse` function that
	transparently handles compressed files.

	Returns iterator over sequence data in file. File is parsed lazily, and so must be kept open.
	The returned iterator is of type :class:`gambit.util.io.ClosingIterator` so it will close the
	file stream automatically when it finishes. It may also be used as a context manager that closes
	the stream on exit. You may also close the stream explicitly using the iterator's ``close``
	method.

	Parameters
	----------
	path
		Path to the file.
	format
		String describing the file format as interpreted by :func:`Bio.SeqIO.parse`.
	compression
		String describing compression method of the file, e.g. ``'gzip'``. ``none`` means no
		compression. Default is to determine compression automatically (can only detect gzip or
		none). See :func:`gambit.util.io.open_compressed`.
	kwargs
		Keyword arguments to :func:`gambit.util.io.open_compressed`.

	Returns
	-------
	gambit.util.io.ClosingIterator
		Iterator yielding :class:`Bio.SeqIO.SeqRecord` instances for each sequence in the file.
	"""

	fobj = open_compressed(path, 'rt', compression, **kwargs)

	try:
		records = SeqIO.parse(fobj, format)
		return ClosingIterator(records, fobj)

	except:
		fobj.close()
		raise
