"""Generic code for working with sequence data.

Note that all code in this package operates on DNA sequences as sequences of
bytes containing ascii-encoded nucleotide codes.

.. data:: NUCLEOTIDES

	``bytes`` corresponding to the four DNA nucleotides. Ascii-encoded upper
	case letters ``ACGT``. Note that the order, while arbitrary, is important
	in this variable as it defines how unique indices are assigned to k-mer
	sequences.
"""

from typing import Union

from Bio.Seq import Seq

from gambit._cython.kmers import revcomp


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
