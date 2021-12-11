"""Test the gambit.seqs module."""

from io import StringIO
from pathlib import Path

import pytest
import numpy as np
from Bio import Seq, SeqIO

from gambit.seq import SequenceFile, revcomp
from gambit.kmers import nkmers, index_to_kmer
import gambit.io.util as ioutil
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

	@pytest.fixture(params=list(ioutil.COMPRESSED_OPENERS), scope='class')
	def compression(self, request):
		"""SequenceFile.compression attribute."""
		return request.param

	@pytest.fixture()
	def info(self, tmpdir, format, compression):
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

	@pytest.fixture
	def info_exists(self, info, seqrecords):
		"""Copy of "info" fixture, but with "seqrecords" written to the file."""

		with info.open('rt') as fobj:
			SeqIO.write(seqrecords, fobj, info.format)

	def test_constructor(self):
		"""Test constructor."""

		info = SequenceFile('foo.fasta', 'fasta')
		assert info == SequenceFile('foo.fasta', 'fasta', None)
		assert info.path == Path('foo.fasta')

	def test_eq(self):
		"""Test equality checking of instances."""
		infos = [
			SequenceFile(p, format, comp)
			for p in ['foo', 'bar']
			for format in ['fasta', 'genbank']
			for comp in [None, 'gzip']
		]

		for i, info1 in enumerate(infos):
			for j, info2 in enumerate(infos):
				if i == j:
					# Try with different instance
					assert info1 == SequenceFile(info1.path, info1.format, info1.compression)
				else:
					assert info1 != info2

	@pytest.mark.parametrize('binary', [False, True])
	def test_open(self, info, file_contents, binary):
		"""Test sequence file is readable and writable."""

		to_write = file_contents.encode() if binary else file_contents

		# Write data to file
		with info.open('wb' if binary else 'wt') as fobj:
			fobj.write(to_write)

		# Read it back and make sure it's the same
		with info.open('rb' if binary else 'rt') as fobj:
			read = fobj.read()

		assert read == to_write

	def test_parse(self, info, seqrecords, file_contents):
		"""Test the parse() method, ensure we get the right records back."""

		# Write pre-formatted contents to file
		with info.open('w') as fobj:
			fobj.write(file_contents)

		# Parse the sequences from it
		parsed = list(info.parse())

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

		info1 = SequenceFile(path, 'fasta')
		assert isinstance(info1, SequenceFile) and info1.path == path

		info2 = SequenceFile(str(path), 'fasta')
		assert isinstance(info2, SequenceFile) and info2.path == path

	def test_absolute(self):
		"""Test the absolute() method."""

		relinfo = SequenceFile('foo/bar.fasta', 'fasta')
		assert not relinfo.path.is_absolute()

		absinfo = relinfo.absolute()
		assert absinfo.path.is_absolute()
		assert absinfo.path == relinfo.path.absolute()

		absinfo2 = absinfo.absolute()
		assert absinfo2 == absinfo

	def test_from_paths(self, format, compression):
		"""Test the from_paths() class method."""

		# List of unique path strings
		paths = ['foo/bar{}.{}'.format(i, format) for i in range(20)]

		infos = SequenceFile.from_paths(paths, format, compression)

		assert len(paths) == len(infos)

		for path, info in zip_strict(paths, infos):
			assert isinstance(info, SequenceFile)
			assert str(info.path) == path
			assert info.format == format
			assert info.compression == compression
