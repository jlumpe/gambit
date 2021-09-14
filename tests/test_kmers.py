"""Tests for gambit.kmers module."""

import pytest
import numpy as np
from Bio.Seq import Seq

from gambit import kmers
from gambit.kmers import KmerSpec
import gambit.io.json as gjson
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

SEQ_TYPES = [str, bytes, bytearray, Seq]

def seq_to_type(seq, type):
	seq = kmers.seq_to_bytes(seq)
	if type is str:
		return seq.decode('ascii')
	return type(seq)


class TestIndices:
	"""Test representation of k-mers by their indices."""

	def test_index_dtype(self):
		"""Test index_dtype() function."""

		# Try k from 0 to 32 (all have dtypes)
		for k in range(33):
			# Check dtype can store the largest index
			top_idx = kmers.nkmers(k) - 1
			assert kmers.index_dtype(k).type(top_idx) == top_idx

		# k > 32 should have no dtype
		assert kmers.index_dtype(33) is None

	def test_nucleotide_order(self):
		"""Check k-mer indices correspond to defined nucleotide order."""

		for i, nuc in enumerate(kmers.NUCLEOTIDES):
			assert kmers.kmer_to_index(bytes([nuc])) == i

	def test_index_conversion(self):
		"""Test converting k-mers to and from their indices."""

		# Test for k in range 0-10
		for k in range(11):

			# Test all indices to max of 1000
			for index in range(min(kmers.nkmers(k), 1000)):

				# Check getting kmer from index
				kmer = kmers.index_to_kmer(index, k)
				assert len(kmer) == k
				assert set(kmers.NUCLEOTIDES).issuperset(kmer)

				# Check conversion back to index
				for T in SEQ_TYPES:
					assert kmers.kmer_to_index(seq_to_type(kmer, T)) == index
					assert kmers.kmer_to_index(seq_to_type(kmer.lower(), T)) == index

					rc = kmers.revcomp(kmer)
					assert kmers.kmer_to_index_rc(seq_to_type(rc, T)) == index
					assert kmers.kmer_to_index_rc(seq_to_type(rc.lower(), T)) == index

		# Check invalid raises error
		with pytest.raises(ValueError):
			kmers.kmer_to_index(b'ATGNC')


class TestKmerSpec:
	"""Test gambit.kmers.KmerSpec."""

	def test_constructor(self):
		# Prefix conversion
		for T in SEQ_TYPES:
			assert KmerSpec(11, seq_to_type('ATGAC', T)).prefix == b'ATGAC'
			assert KmerSpec(11, seq_to_type('atgac', T)).prefix == b'ATGAC'

		# Invalid k
		with pytest.raises(ValueError):
			KmerSpec(0, 'ATGAC')

		# Invalid prefix
		for prefix in [b'ATGAX', 'ATGAX']:
			with pytest.raises(ValueError):
				KmerSpec(11, prefix)

	def test_attributes(self):
		"""Test basic attributes."""

		# Try k from 1 to 32 (all have dtypes)
		for k in range(1, 33):

			spec = KmerSpec(k, 'ATGAC')

			# Check length attributes
			assert spec.prefix_len == len(spec.prefix)
			assert spec.prefix_len + spec.k == spec.total_len

		# Check prefix is bytes
		assert isinstance(KmerSpec(11, 'ATGAC').prefix, bytes)

	def test_eq(self):
		"""Test equality testing."""

		kspec = KmerSpec(11, 'ATGAC')
		assert kspec == KmerSpec(11, 'ATGAC')
		assert kspec != KmerSpec(11, 'ATGAA')
		assert kspec != KmerSpec(12, 'ATGAC')

	def test_pickle(self):

		import pickle

		kspec = KmerSpec(11, 'ATGAC')

		assert kspec == pickle.loads(pickle.dumps(kspec))

	def test_json(self):
		"""Test conversion to/from JSON."""

		kspec = KmerSpec(11, 'ATGAC')
		data = gjson.to_json(kspec)

		assert data == dict(
			k=kspec.k,
			prefix=kspec.prefix.decode('ascii'),
		)

		assert gjson.from_json(data, KmerSpec) == kspec


def test_dense_sparse_conversion():
	"""Test conversion between dense and sparse representations of k-mer coordinates."""

	for k in range(1, 10):

		kspec = KmerSpec(k, 'ATGAC')

		# Create vector with every 3rd k-mer
		vec = np.zeros(kspec.nkmers, dtype=bool)
		vec[np.arange(vec.size) % 3 == 0] = True

		# Convert to coords
		coords = kmers.dense_to_sparse(vec)

		# Check coords
		assert len(coords) == vec.sum()
		for index in coords:
			assert vec[index]

		# Check coords ascending
		assert np.all(np.diff(coords) > 0)

		# Check converting back
		assert np.array_equal(vec, kmers.sparse_to_dense(kspec, coords))


def check_revcomp(seq, rc):
	"""Assert the reverse complement of a sequence is correct."""
	l = len(seq)
	for i in range(l):
		assert rc[l - i - 1] == NUC_COMPLEMENTS.get(seq[i], seq[i])


def test_revcomp():
	"""Test gambit._cython.kmers.revcomp."""

	# Check empty
	assert kmers.revcomp(b'') == b''

	# Check one-nucleotide values
	for nuc1, nuc2 in NUC_COMPLEMENTS.items():
		b1, b2 = [bytes([n]) for n in [nuc1, nuc2]]
		assert kmers.revcomp(b1) == b2
		assert kmers.revcomp(b1.lower()) == b2.lower()

	# Check single invalid code
	assert kmers.revcomp(b'N') == b'N'
	assert kmers.revcomp(b'n') == b'n'

	# Check all 6-mers
	k = 6
	for i in range(kmers.nkmers(k)):
		kmer = kmers.index_to_kmer(i, k)

		rc = kmers.revcomp(kmer)

		check_revcomp(rc, kmer)
		check_revcomp(rc.lower(), kmer.lower())

		assert kmers.revcomp(rc) == kmer
		assert kmers.revcomp(rc.lower()) == kmer.lower()

	# Check longer seqs with invalid nucleotides
	seq = bytearray(b'ATGCatgc')

	for i in range(len(seq)):

		array = bytearray(seq)
		array[i] = ord(b'N')
		seq2 = bytes(array)

		rc = kmers.revcomp(seq2)

		check_revcomp(rc, seq2)
		assert kmers.revcomp(rc) == seq2


class TestKmerSpecConversion:
	"""Test converting signatures from one KmerSpec to another."""

	def test_can_convert(self):
		from_kspec = KmerSpec(11, 'ATGAC')

		compatible = [
			KmerSpec(11, 'ATGAC'),
			KmerSpec(8, 'ATGAC'),
			KmerSpec(10, 'ATGACA'),
			KmerSpec(8, 'ATGACA'),
		]

		for to_kspec in compatible:
			assert kmers.can_convert(from_kspec, to_kspec)
			kmers.check_can_convert(from_kspec, to_kspec)

		incompatible = [
			KmerSpec(11, 'CAGTA'),
			KmerSpec(12, 'ATGAC'),
			KmerSpec(11, 'ATGA'),
			KmerSpec(11, 'ATGACT'),
		]

		for to_kspec in incompatible:
			assert not kmers.can_convert(from_kspec, to_kspec)
			with pytest.raises(ValueError):
				kmers.check_can_convert(from_kspec, to_kspec)

	@pytest.fixture(scope='class')
	def seqs(self):
		np.random.seed(0)
		return [random_seq(100_000) for _ in range(100)]

	@pytest.mark.parametrize('to_kspec', [
		KmerSpec(10, 'ATGAC'),   # Reduce k
		KmerSpec(8, 'ATGAC'),    # Reduce k
		KmerSpec(9, 'ATGACGT'),  # Extend prefix
		KmerSpec(7, 'ATGACGT'),  # Extend prefix and reduce k further
	])
	def test_convert(self, seqs, to_kspec):
		from gambit.search import calc_signature

		from_kspec = KmerSpec(11, 'ATGAC')

		for seq in seqs:
			from_vec = calc_signature(from_kspec, seq, sparse=False)
			from_sig = kmers.dense_to_sparse(from_vec)

			to_vec = kmers.convert_dense(from_kspec, to_kspec, from_vec)
			to_sig = kmers.convert_sparse(from_kspec, to_kspec, from_sig)

			found_vec = calc_signature(to_kspec, seq, sparse=False)

			assert np.array_equal(to_vec, found_vec)
			assert np.array_equal(to_sig, kmers.dense_to_sparse(found_vec))
