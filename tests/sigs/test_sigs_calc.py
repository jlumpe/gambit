"""Tests for gambit.search module."""

from io import StringIO

import pytest
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq

from gambit.sigs.calc import calc_signature, calc_file_signature, calc_file_signatures
from gambit.kmers import KmerSpec, index_to_kmer
from gambit.seq import SEQ_TYPES, revcomp
from gambit.test import fill_bytearray, make_kmer_seq, make_kmer_seqs, check_progress, convert_seq
from gambit.io.seq import SequenceFile
import gambit.io.util as ioutil
from gambit.sigs import sigarray_eq


KSPEC = KmerSpec(11, 'AGTAC')


class TestCalcSignature:
	"""Test the calc_signature() function."""

	@pytest.mark.parametrize('seq_type', SEQ_TYPES)
	def test_single(self, seq_type):
		"""Test with single signature as argument."""

		np.random.seed(0)
		seq_bytes, expected = make_kmer_seq(KSPEC, 100000, 50, 10)
		seq = convert_seq(seq_bytes, seq_type)

		# Test normal
		result = calc_signature(KSPEC, seq)
		assert np.array_equal(result, expected)

		# Test reverse complement
		rcseq = convert_seq(revcomp(seq_bytes), seq_type)
		result = calc_signature(KSPEC, rcseq)
		assert np.array_equal(result, expected)

		# Test lower case
		result = calc_signature(KSPEC, seq.lower())
		assert np.array_equal(result, expected)

	@pytest.mark.parametrize('seq_type', SEQ_TYPES)
	def test_multiple(self, seq_type):
		"""Test with list signatures as argument."""

		np.random.seed(0)
		seqs, sig = make_kmer_seqs(KSPEC, 10, 10000, 50, 10)
		seqs2 = [convert_seq(seq, seq_type) for seq in seqs]

		result = calc_signature(KSPEC, seqs2)
		assert np.array_equal(result, sig)

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
		found = calc_signature(kspec, seq)

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
			revcomp(b'GGATTATATAG'),
			revcomp(b'TCCGATAGGAT'),
		}

		for s in [seq, revcomp(seq)]:
			sig = calc_signature(kspec, s)
			found = [index_to_kmer(idx, kspec.k) for idx in sig]

			assert len(found) == len(expected)
			assert all(kmer in expected for kmer in found)


class TestCalcFileSignatures:

	@pytest.fixture(scope='class')
	def record_sets(self):

		items = []

		np.random.seed(0)
		for i in range(5):
			seqs, sig = make_kmer_seqs(KSPEC, 10, 10000, 50, 10)

			# Create the BioPython sequence record object
			records = [SeqIO.SeqRecord(
				seq=Seq(seq.decode('ascii')),
				id='SEQ{}'.format(i + 1),
				description='sequence {}'.format(i + 1),
			) for seq in seqs]

			items.append((records, sig))

		return items

	@pytest.fixture(scope='class', params=['fasta'])
	def format(self, request):
		return request.param

	@pytest.fixture(scope='class', params=list(ioutil.COMPRESSED_OPENERS))
	def compression(self, request):
		return request.param

	@pytest.fixture()
	def files(self, record_sets, tmp_path, format, compression):

		files = []

		for i, (records, sig) in enumerate(record_sets):
			file = SequenceFile(tmp_path / f'{i + 1}.fasta', format, compression)

			with file.open('w') as f:
				SeqIO.write(records, f, format)

			files.append(file)

		return files

	def test_calc_file_signature(self, record_sets, files):
		"""Test the calc_file_signature function."""

		for file, (records, sig) in zip(files, record_sets):
			result = calc_file_signature(KSPEC, file)
			assert np.array_equal(result, sig)

	@pytest.mark.parametrize('concurrency', [None, 'threads', 'processes'])
	def test_calc_file_signatures(self, record_sets, files, concurrency):
		"""Test the calc_file_signatures function."""
		sigs = [sig for records, sig in record_sets]

		with check_progress(total=len(files)) as pconf:
			sigs2 = calc_file_signatures(KSPEC, files, progress=pconf, concurrency=concurrency)

		assert sigarray_eq(sigs, sigs2)
