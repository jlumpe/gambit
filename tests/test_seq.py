"""Test the gambit.seqs module."""

from gambit.seq import revcomp
from gambit.kmers import nkmers, index_to_kmer


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
