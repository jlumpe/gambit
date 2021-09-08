"""Tests for gambit.search module."""

from io import StringIO

import pytest
import numpy as np
from Bio import Seq, SeqIO

from gambit.search import find_kmers, find_kmers_parse, find_kmers_in_file, find_kmers_in_files
from gambit.kmers import KmerSpec, reverse_complement, dense_to_sparse, sparse_to_dense, index_to_kmer
from gambit.test import fill_bytearray, make_kmer_seq, check_progress

from gambit.io.seq import SequenceFile
import gambit.io.util as ioutil
from gambit.signatures import sigarray_eq


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
	vec = np.zeros(kspec.nkmers, dtype=bool)

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
	"""Test the find_kmers() function."""

	@pytest.mark.parametrize('sparse', [True, False])
	def test_basic(self, sparse):
		"""Test general k-mer finding."""

		kspec = KmerSpec(11, 'ATGAC')

		np.random.seed(0)
		seq, signature = make_kmer_seq(kspec, 100000, kmer_interval=50, n_interval=10)
		expected = signature if sparse else sparse_to_dense(kspec, signature)

		# Test normal
		result = find_kmers(kspec, seq, sparse=sparse)
		assert np.array_equal(result, expected)

		# Test reverse complement
		result = find_kmers(kspec, reverse_complement(seq), sparse=sparse)
		assert np.array_equal(result, expected)

		# Test lower case
		result = find_kmers(kspec, seq.lower(), sparse=sparse)
		assert np.array_equal(result, expected)

		# Test string argument
		result = find_kmers(kspec, seq.decode('ascii'), sparse=sparse)
		assert np.array_equal(result, expected)

	def test_bounds(self):
		"""Test k-mer finding at beginning and end of sequence to catch errors with search bounds."""

		# Sequence of all ATN's
		seqlen = 100000
		seq_array = fill_bytearray(b'ATN', seqlen)

		# Choose prefix with nucleotides not found in sequence "background"
		kspec = KmerSpec(11, b'CCGGG')

		# Add at beginning
		seq_array[0:kspec.prefix_len] = kspec.prefix
		seq_array[kspec.prefix_len:kspec.total_len] = index_to_kmer(0, kspec.k)

		# Add at end
		seq_array[-kspec.total_len:-kspec.k] = kspec.prefix
		seq_array[-kspec.k:] = index_to_kmer(1, kspec.k)

		seq = bytes(seq_array)
		found = find_kmers(kspec, seq)

		assert np.array_equal(found, [0, 1])

	def test_overlapping(self):
		"""Test k-mer finding when k-mers overlap with each other.

		The test sequence is manually designed to have a variety of overlapping
		forwards and backwards matches
		"""

		kspec = KmerSpec(11, b'GCCGG')

		seq = b'ATATGCCGGCCGGATTATATAGCCGGCATTACATCCGATAGGATCCGGCAATAA'
		#      |    |>>>>...........
		#      |        |>>>>........... (forward match which overlaps prefix)
		#      |                     |>>>>........... (another overlapping forward match)
		#      |....<<<<| (backward match for prefix, but too close to end)
		#      |           ...........<<<<|
		#      |                                 ...........<<<<|

		expected = {
			b'CCGGATTATAT',
			b'ATTATATAGCC',
			b'CATTACATCCG',
			reverse_complement(b'GGATTATATAG'),
			reverse_complement(b'TCCGATAGGAT'),
		}

		for s in [seq, reverse_complement(seq)]:
			sig = find_kmers(kspec, s)
			found = [index_to_kmer(idx, kspec.k) for idx in sig]

			assert len(found) == len(expected)
			assert all(kmer in expected for kmer in found)


class TestFindKmersInFile:
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

	@pytest.mark.parametrize('concurrency', [None, 'threads', 'processes'])
	def test_find_kmers_in_files(self, seq_data, files, concurrency):
		"""Test the find_kmers_in_files function."""
		seqs, sigs = seq_data

		with check_progress(total=len(files)) as pconf:
			sigs2 = find_kmers_in_files(self.KSPEC, files, progress=pconf, concurrency=concurrency)

		assert sigarray_eq(sigs, sigs2)
