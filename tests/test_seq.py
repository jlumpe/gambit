"""Test the gambit.seqs module."""

from io import StringIO
from pathlib import Path
import os

import pytest
import numpy as np
from Bio import Seq, SeqIO

from gambit.seq import revcomp, parse_seqs
from gambit.kmers import nkmers, index_to_kmer
from gambit.util.misc import zip_strict
from gambit.util.io import open_compressed

from .common import random_seq


# Complements to nucleotide ASCII codes
NUC_COMPLEMENTS = {
	65: 84,
	84: 65,
	71: 67,
	67: 71,
	97: 116,
	116: 97,
	103: 99,
	99: 103,
}


def check_revcomp(seq, rc):
	"""Assert the reverse complement of a sequence is correct."""
	l = len(seq)
	for i in range(l):
		assert rc[l - i - 1] == NUC_COMPLEMENTS.get(seq[i], seq[i])


def test_revcomp():
	"""Test gambit._cython.revcomp."""

	# Check empty
	assert revcomp(b'') == b''

	# Check one-nucleotide values
	for nuc1, nuc2 in NUC_COMPLEMENTS.items():
		b1, b2 = [bytes([n]) for n in [nuc1, nuc2]]
		assert revcomp(b1) == b2
		assert revcomp(b1.lower()) == b2.lower()

	# Check single invalid code
	assert revcomp(b'N') == b'N'
	assert revcomp(b'n') == b'n'

	# Check all 6-mers
	k = 6
	for i in range(nkmers(k)):
		kmer = index_to_kmer(i, k)

		rc = revcomp(kmer)

		check_revcomp(rc, kmer)
		check_revcomp(rc.lower(), kmer.lower())

		assert revcomp(rc) == kmer
		assert revcomp(rc.lower()) == kmer.lower()

	# Check longer seqs with invalid nucleotides
	seq = bytearray(b'ATGCatgc')

	for i in range(len(seq)):

		array = bytearray(seq)
		array[i] = ord(b'N')
		seq2 = bytes(array)

		rc = revcomp(seq2)

		check_revcomp(rc, seq2)
		assert revcomp(rc) == seq2


def test_seq_to_bytes():
	"""Test seq_to_bytes() function."""
	# TODO


def test_validate_dna_seq_bytes():
	"""Test validate_dna_seq_bytes() function."""
	# TODO


@pytest.fixture(scope='module')
def seqrecords():
	"""Random SeqRecord instances."""

	records = []

	np.random.seed(0)

	for i in range(20):
		seq = Seq.Seq(random_seq(1000).decode('ascii'))
		id_ = f'seq{i + 1}'
		descr = f'{id_} Test sequence {i + 1}'
		records.append(SeqIO.SeqRecord(seq, id=id_, description=descr))

	return records


@pytest.mark.parametrize('compression', ['none', 'gzip'])
@pytest.mark.parametrize('auto', [False, True])
def test_parse_seqs(tmp_path: Path, seqrecords: list[SeqIO.SeqRecord], compression: str, auto: bool):
	"""Test the parse_seqs() function."""

	# Write FASTA file
	file = tmp_path / ('test.fa' + ('.gz' if compression == 'gzip' else ''))
	with open_compressed(file, 'wt', compression) as fh:
		SeqIO.write(seqrecords, fh, 'fasta')

	# Parse
	with parse_seqs(file, 'fasta', compression='auto' if auto else compression) as parsed:
		records2 = list(parsed)

		# Check ClosingIterator closes the underlying file object when last record is read
		assert parsed.fobj.closed

	# Check parsed records are correct
	for record, record2 in zip_strict(seqrecords, records2):
		assert record2.seq == record.seq
		assert record2.id == record.id
		assert record2.description == record.description
