"""Tests for gambit.search module."""

from io import StringIO

import pytest
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq

from gambit.search import find_kmers, calc_signature, calc_signature_parse, calc_file_signature, \
	calc_file_signatures
from gambit.kmers import KmerSpec, reverse_complement, dense_to_sparse, sparse_to_dense, \
	index_to_kmer, kmer_to_index, NUCLEOTIDES
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
			seq=Seq(seq.decode('ascii')),
			id='SEQ{}'.format(i + 1),
			description='sequence {}'.format(i + 1),
		))

	return records, vec


def convert_seq(seq, t: type):
	if isinstance(seq, t):
		return seq

	if t is bytes:
		return seq.encode('ascii')

	if isinstance(seq, bytes):
		seq = seq.decode('ascii')

	if t is str:
		return str(seq)
	elif t is Seq:
		return Seq(seq)
	else:
		assert 0


@pytest.fixture(params=[bytes, str, Seq])
def seq_type(request):
	return request.param


def test_find_kmers(seq_type):
	"""Test the find_kmers() function and KmerMatch class."""

	kspec = KmerSpec(11, 'ATGAC')

	np.random.seed(0)
	seq, sig = make_kmer_seq(kspec, 100000, kmer_interval=50, n_interval=10)
	seq = convert_seq(seq, seq_type)

	found = []

	for match in find_kmers(kspec, seq):
		assert match.kmerspec is kspec
		assert match.seq is seq

		kmer_indices = match.kmer_indices()
		full_indices = match.full_indices()

		assert kmer_indices.stop - kmer_indices.start == kspec.k
		assert full_indices.stop - full_indices.start == kspec.total_len
		if match.reverse:
			assert full_indices.start == kmer_indices.start
		else:
			assert full_indices.stop == kmer_indices.stop

		matched = convert_seq(seq[kmer_indices], bytes)
		if match.reverse:
			matched = reverse_complement(matched)

		matched_full = convert_seq(seq[full_indices], bytes)
		if match.reverse:
			matched_full = reverse_complement(matched_full)

		assert matched_full == kspec.prefix + matched
		assert match.kmer() == matched

		try:
			index = kmer_to_index(matched)
		except ValueError:
			assert any(c not in NUCLEOTIDES for c in matched.upper())
			continue

		assert match.kmer_index() == index
		found.append(index)

	assert np.array_equal(sorted(found), sig)


@pytest.mark.parametrize('sparse', [True, False])
def test_calc_signature(sparse, seq_type):
	"""Test the calc_signature() function."""

	kspec = KmerSpec(11, 'ATGAC')

	np.random.seed(0)
	seq_bytes, signature = make_kmer_seq(kspec, 100000, kmer_interval=50, n_interval=10)
	seq = convert_seq(seq_bytes, seq_type)
	expected = signature if sparse else sparse_to_dense(kspec, signature)

	# Test normal
	result = calc_signature(kspec, seq, sparse=sparse)
	assert np.array_equal(result, expected)

	# Test reverse complement
	rcseq = convert_seq(reverse_complement(seq_bytes), seq_type)
	result = calc_signature(kspec, rcseq, sparse=sparse)
	assert np.array_equal(result, expected)

	# Test lower case
	result = calc_signature(kspec, seq.lower(), sparse=sparse)
	assert np.array_equal(result, expected)


def test_bounds():
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
	found = calc_signature(kspec, seq)

	assert np.array_equal(found, [0, 1])


def test_overlapping():
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
		sig = calc_signature(kspec, s)
		found = [index_to_kmer(idx, kspec.k) for idx in sig]

		assert len(found) == len(expected)
		assert all(kmer in expected for kmer in found)


class TestCalcFileSignatures:
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
	def test_calc_signature_parse(self, seq_data, format, sparse):
		"""Test the calc_signature_parse function."""

		for records, sig in zip(*seq_data):
			# Parse from buffer
			buf = StringIO()
			SeqIO.write(records, buf, format)
			buf.seek(0)

			result = calc_signature_parse(self.KSPEC, buf, 'fasta', sparse=sparse)

			if sparse:
				assert np.array_equal(result, sig)
			else:
				assert np.array_equal(dense_to_sparse(result), sig)

	@pytest.mark.parametrize('sparse', [False, True])
	def test_calc_file_signature(self, seq_data, files, sparse):
		"""Test the calc_file_signature function."""

		seqs, sigs = seq_data

		for file, sig in zip(files, sigs):
			result = calc_file_signature(self.KSPEC, file, sparse=sparse)

			if sparse:
				assert np.array_equal(result, sig)
			else:
				assert np.array_equal(dense_to_sparse(result), sig)

	@pytest.mark.parametrize('concurrency', [None, 'threads', 'processes'])
	def test_calc_file_signatures(self, seq_data, files, concurrency):
		"""Test the calc_file_signatures function."""
		seqs, sigs = seq_data

		with check_progress(total=len(files)) as pconf:
			sigs2 = calc_file_signatures(self.KSPEC, files, progress=pconf, concurrency=concurrency)

		assert sigarray_eq(sigs, sigs2)
