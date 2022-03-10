"""Test the gambit.seqs module."""

from io import StringIO
from pathlib import Path
import os

import pytest
import numpy as np
from Bio import Seq, SeqIO

from gambit.seq import SequenceFile, revcomp
from gambit.kmers import nkmers, index_to_kmer
from gambit.util.misc import zip_strict
from gambit.test import random_seq


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


class TestSequenceFile:
	"""Test the SequenceFile class."""

	@pytest.fixture(params=['fasta'], scope='class')
	def format(self, request):
		"""SequenceFile.format attribute."""
		return request.param

	@pytest.fixture(params=[None, 'gzip'], scope='class')
	def compression(self, request):
		"""SequenceFile.compression attribute."""
		return request.param

	@pytest.fixture()
	def seqfile(self, tmpdir, format, compression):
		"""A SequenceFile instance pointing to a file in a test temporary directory.

		File does not yet exist.
		"""
		path = tmpdir.join('test.' + format).strpath
		return SequenceFile(path, format, compression)

	@pytest.fixture(scope='class')
	def seqrecords(self):
		"""A collection of random Bio.SeqIO.SeqRecord's."""
		np.random.seed(0)
		records = []

		for i in range(20):
			seq = Seq.Seq(random_seq(1000).decode('ascii'))
			id_ = 'seq{}'.format(i + 1)
			descr = 'Test sequence {}'.format(i + 1)
			records.append(SeqIO.SeqRecord(seq, id=id_, description=descr))

		return tuple(records)

	@pytest.fixture
	def file_contents(self, format, seqrecords):
		"""String contents of a file containing the sequence records."""
		buf = StringIO()
		SeqIO.write(seqrecords, buf, format)
		return buf.getvalue()

	def test_constructor(self):
		"""Test constructor."""

		seqfile = SequenceFile('foo.fasta', 'fasta')
		assert seqfile == SequenceFile('foo.fasta', 'fasta', None)
		assert seqfile.path == Path('foo.fasta')

	def test_eq(self):
		"""Test equality checking of instances."""
		seqfiles = [
			SequenceFile(p, format, comp)
			for p in ['foo', 'bar']
			for format in ['fasta', 'genbank']
			for comp in [None, 'gzip']
		]

		for i, seqfile1 in enumerate(seqfiles):
			for j, seqfile2 in enumerate(seqfiles):
				if i == j:
					# Try with different instance
					assert seqfile1 == SequenceFile(seqfile1.path, seqfile1.format, seqfile1.compression)
				else:
					assert seqfile1 != seqfile2

	def test_special_methods(self, seqfile):
		assert str(seqfile) == str(seqfile.path)
		assert os.fspath(seqfile) == str(seqfile)

		# Check os.PathLike interface
		text = 'foo'
		with open(seqfile, 'w') as f:
			f.write(text)
		with open(seqfile, 'r') as f:
			read = f.read()
		assert read == text

	@pytest.mark.parametrize('binary', [False, True])
	def test_open(self, seqfile, file_contents, binary):
		"""Test sequence file is readable and writable."""

		to_write = file_contents.encode() if binary else file_contents

		# Write data to file
		with seqfile.open('wb' if binary else 'wt') as fobj:
			fobj.write(to_write)

		# Read it back and make sure it's the same
		with seqfile.open('rb' if binary else 'rt') as fobj:
			read = fobj.read()

		assert read == to_write

	def test_parse(self, seqfile, seqrecords, file_contents):
		"""Test the parse() method, ensure we get the right records back."""

		# Write pre-formatted contents to file
		with seqfile.open('wt') as fobj:
			fobj.write(file_contents)

		# Parse the sequences from it
		parsed = list(seqfile.parse())

		# Check they match
		assert len(parsed) == len(seqrecords)

		for parsed_req, orig_req in zip_strict(parsed, seqrecords):
			assert isinstance(parsed_req, SeqIO.SeqRecord)
			assert parsed_req.seq == orig_req.seq
			assert parsed_req.id == orig_req.id

			# This is something stupid BioPython does - when writing a SeqRecord
			# as FASTA it writes the .id attributed followed by a space and then
			# the .description attribute on the description line. When reading,
			# the entire line is used as the description attribute and so
			# includes the ID
			assert parsed_req.description == orig_req.id + ' ' + orig_req.description

	def test_path_arg(self):
		"""Test the "path" argument to the constructor."""

		path = Path('foo/bar.fasta')

		seqfile1 = SequenceFile(path, 'fasta')
		assert isinstance(seqfile1, SequenceFile) and seqfile1.path == path

		seqfile2 = SequenceFile(str(path), 'fasta')
		assert isinstance(seqfile2, SequenceFile) and seqfile2.path == path

	def test_absolute(self):
		"""Test the absolute() method."""

		relseqfile = SequenceFile('foo/bar.fasta', 'fasta')
		assert not relseqfile.path.is_absolute()

		absseqfile = relseqfile.absolute()
		assert absseqfile.path.is_absolute()
		assert absseqfile.path == relseqfile.path.absolute()

		absseqfile2 = absseqfile.absolute()
		assert absseqfile2 == absseqfile

	def test_from_paths(self, format, compression):
		"""Test the from_paths() class method."""

		# List of unique path strings
		paths = ['foo/bar{}.{}'.format(i, format) for i in range(20)]

		seqfiles = SequenceFile.from_paths(paths, format, compression)

		assert len(paths) == len(seqfiles)

		for path, seqfile in zip_strict(paths, seqfiles):
			assert isinstance(seqfile, SequenceFile)
			assert str(seqfile.path) == path
			assert seqfile.format == format
			assert seqfile.compression == compression
