"""Test gambit.io.seq."""

from io import StringIO
from pathlib import Path

import pytest
import numpy as np
from Bio import Seq, SeqIO

from gambit.io.seq import SequenceFile, find_kmers_parse, find_kmers_in_file, find_kmers_in_files
import gambit.io.util as ioutil
from gambit.kmers import KmerSpec, dense_to_sparse, sparse_to_dense
from gambit.signatures import sigarray_eq
from gambit.util.misc import zip_strict
from gambit.test import make_kmer_seq, random_seq


def create_sequence_records(kspec, n, seq_len=10000):
	"""Create a set of random DNA sequences with known combined k-mer signature.

	Parameters
	----------
	kspec : KmerSpec
	n : int
		Number of sequences to create.
	seq_len : int
		Length of sequences to create.

	Returns
	-------
	tuple
		(records, kmer_vec) tuple.
	"""
	records = []
	vec = np.zeros(4 ** kspec.k, dtype=bool)

	for i in range(n):
		seq, sig = make_kmer_seq(kspec, seq_len, kmer_interval=50, n_interval=10)

		# Combine vectors of all sequences
		vec |= sparse_to_dense(kspec, sig)

		# Convert every other sequence to lower case, just to switch things up...
		if i % 2:
			seq = seq.lower()

		# Create the BioPython sequence record object
		records.append(SeqIO.SeqRecord(
			seq=Seq.Seq(seq.decode('ascii')),
			id='SEQ{}'.format(i + 1),
			description='sequence {}'.format(i + 1),
		))

	return records, vec


class TestFindKmers:
	KSPEC = KmerSpec(11, 'AGTAC')

	@pytest.fixture(scope='class')
	def seq_data(self):
		n = 5

		seqs = []
		sigs = []

		# Create files
		np.random.seed(0)
		for i in range(n):
			records, vec = create_sequence_records(self.KSPEC, 10)
			seqs.append(records)
			sigs.append(dense_to_sparse(vec))

		return seqs, sigs

	@pytest.fixture(scope='class', params=['fasta'])
	def format(self, request):
		return request.param

	@pytest.fixture(scope='class', params=list(ioutil.COMPRESSED_OPENERS))
	def compression(self, request):
		return request.param

	@pytest.fixture()
	def files(self, seq_data, tmp_path, format, compression):
		seqs, sigs = seq_data

		files = []

		for i, records in enumerate(seqs):
			file = SequenceFile(tmp_path / f'{i + 1}.fasta', format, compression)

			with file.open('w') as f:
				SeqIO.write(records, f, format)

			files.append(file)

		return files

	@pytest.mark.parametrize('sparse', [False, True])
	def test_find_kmers_parse(self, seq_data, format, sparse):
		"""Test the find_kmers_parse function."""

		for records, sig in zip(*seq_data):
			# Parse from buffer
			buf = StringIO()
			SeqIO.write(records, buf, format)
			buf.seek(0)

			result = find_kmers_parse(self.KSPEC, buf, 'fasta', sparse=sparse)

			if sparse:
				assert np.array_equal(result, sig)
			else:
				assert np.array_equal(dense_to_sparse(result), sig)

	@pytest.mark.parametrize('sparse', [False, True])
	def test_find_kmers_in_file(self, seq_data, files, sparse):
		"""Test the find_kmers_in_file function."""

		seqs, sigs = seq_data

		for file, sig in zip(files, sigs):
			result = find_kmers_in_file(self.KSPEC, file, sparse=sparse)

			if sparse:
				assert np.array_equal(result, sig)
			else:
				assert np.array_equal(dense_to_sparse(result), sig)

	def test_find_kmers_in_files(self, seq_data, files):
		"""Test the find_kmers_in_files function.

		TODO - test progress monitoring
		"""
		seqs, sigs = seq_data

		sigs2 = find_kmers_in_files(self.KSPEC, files)
		assert sigarray_eq(sigs, sigs2)


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
