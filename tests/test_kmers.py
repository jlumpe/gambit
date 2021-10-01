"""Tests for gambit.kmers module."""

import pytest
import numpy as np

from gambit.seq import SEQ_TYPES, NUCLEOTIDES, revcomp
from gambit import kmers
from gambit.kmers import KmerSpec
import gambit.io.json as gjson
from gambit.test import convert_seq, make_kmer_seq


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

		for i, nuc in enumerate(NUCLEOTIDES):
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
				assert set(NUCLEOTIDES).issuperset(kmer)

				# Check conversion back to index
				for T in SEQ_TYPES:
					assert kmers.kmer_to_index(convert_seq(kmer, T)) == index
					assert kmers.kmer_to_index(convert_seq(kmer.lower(), T)) == index

					rc = revcomp(kmer)
					assert kmers.kmer_to_index_rc(convert_seq(rc, T)) == index
					assert kmers.kmer_to_index_rc(convert_seq(rc.lower(), T)) == index

		# Check invalid raises error
		with pytest.raises(ValueError):
			kmers.kmer_to_index(b'ATGNC')


class TestKmerSpec:
	"""Test gambit.kmers.KmerSpec."""

	def test_constructor(self):
		# Prefix conversion
		for T in SEQ_TYPES:
			assert KmerSpec(11, convert_seq('ATGAC', T)).prefix == b'ATGAC'
			assert KmerSpec(11, convert_seq('atgac', T)).prefix == b'ATGAC'

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
		assert hash(kspec) == hash(KmerSpec(11, 'ATGAC'))
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


@pytest.mark.parametrize('lower', [False, True])
@pytest.mark.parametrize('seq_type', SEQ_TYPES)
def test_find_kmers(seq_type, lower):
	"""Test the find_kmers() function and KmerMatch class."""

	kspec = KmerSpec(11, 'ATGAC')

	np.random.seed(0)
	seq, sig = make_kmer_seq(kspec, 100000, kmer_interval=50, n_interval=10)

	seq = convert_seq(seq, seq_type)
	if lower:
		seq = seq.lower()

	found = []

	for match in kmers.find_kmers(kspec, seq):
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

		matched = convert_seq(seq[kmer_indices], bytes).upper()
		if match.reverse:
			matched = revcomp(matched)

		matched_full = convert_seq(seq[full_indices], bytes).upper()
		if match.reverse:
			matched_full = revcomp(matched_full)

		assert matched_full == kspec.prefix + matched
		assert match.kmer().upper() == matched

		try:
			index = kmers.kmer_to_index(matched)
		except ValueError:
			assert any(c not in NUCLEOTIDES for c in matched.upper())
			continue

		assert match.kmer_index() == index
		found.append(index)

	assert np.array_equal(sorted(found), sig)
