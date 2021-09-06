"""Test gambit.test module."""

import pytest
import numpy as np

from gambit.test import make_signatures, random_seq, fill_bytearray, make_kmer_seq, check_progress
from gambit.kmers import KmerSpec, reverse_complement, kmer_to_index, dense_to_sparse, nkmers
from gambit.util.progress import get_progress


@pytest.mark.parametrize('k', [4, 6, 8])
@pytest.mark.parametrize('n', [10, 100])
@pytest.mark.parametrize('dtype', [np.dtype('u8'), np.dtype('u4')])
def test_make_signatures(k, n, dtype):
	np.random.seed(0)
	sigs = make_signatures(k, n, dtype)
	assert len(sigs) == n

	for i, sig in enumerate(sigs):
		assert sig.dtype == dtype
		assert np.all(np.diff(sig) > 0)  # sorted
		assert np.all(sig < nkmers(k))  # in expected range

		# Pairwise distinct
		for j in range(i+1, n):
			assert not np.array_equal(sig, sigs[j])


@pytest.mark.parametrize('n', [100, 1000])
@pytest.mark.parametrize('chars', ['ACGT', 'XYZ'])
def test_random_seq(n, chars):
	np.random.seed(0)
	seq = random_seq(n, chars)
	assert isinstance(seq, bytes)
	assert len(seq) == n
	assert all(chr(c) in chars for c in seq)


@pytest.mark.parametrize('pattern', [b'N', b'ABC'])
@pytest.mark.parametrize('n', [100, 1000])
def test_fill_bytearray(pattern, n):
	arr = fill_bytearray(pattern, n)
	assert isinstance(arr, bytearray)
	assert len(arr) == n

	for i, b in enumerate(arr):
		assert b == pattern[i % len(pattern)]


@pytest.mark.parametrize('kspec', [KmerSpec(11, 'ATGAC'), KmerSpec(8, 'ATG')])
@pytest.mark.parametrize('seqlen', [1_000, 10_000])
@pytest.mark.parametrize('kmer_interval', [20, 50])
@pytest.mark.parametrize('n_interval', [None, 5])
def test_make_kmer_seq(kspec, seqlen, kmer_interval, n_interval):
	np.random.seed(0)
	seq, sig = make_kmer_seq(kspec, seqlen, kmer_interval, n_interval)
	assert len(seq) == seqlen

	vec = np.zeros(kspec.nkmers, dtype=bool)

	for i, p in enumerate(range(0, seqlen - kspec.total_len, kmer_interval)):
		match = seq[p:(p + kspec.total_len)]

		fwd = match.startswith(kspec.prefix)
		rev = match.endswith(reverse_complement(kspec.prefix))
		assert fwd ^ rev

		if rev:
			match = reverse_complement(match)

		kmer = match[kspec.prefix_len:]
		if b'N' not in kmer:
			vec[kmer_to_index(kmer)] = True

	assert np.array_equal(sig, dense_to_sparse(vec))


def test_check_progress():
	"""Test the check_progress function."""

	with check_progress() as pconf:
		with get_progress(pconf, 100) as meter:
			meter.moveto(100)

	with check_progress(total=100) as pconf:
		with get_progress(pconf, 100) as meter:
			meter.moveto(100)

	with check_progress(check_closed=False) as pconf:
		meter = get_progress(pconf, 100)
		meter.moveto(100)

	# Not completed
	with pytest.raises(AssertionError):
		with check_progress() as pconf:
			with get_progress(pconf, 100) as meter:
				meter.moveto(99)

	# Wrong total
	with pytest.raises(AssertionError):
		with check_progress(total=100) as pconf:
			with get_progress(pconf, 10) as meter:
				meter.moveto(10)

	# Not closed
	with pytest.raises(AssertionError):
		with check_progress() as pconf:
			meter = get_progress(pconf, 100)
			meter.moveto(100)

	# Not instantiated
	with pytest.raises(AssertionError):
		with check_progress():
			pass

	# Instantiated multiple times
	with pytest.raises(AssertionError):
		with check_progress() as pconf:
			with get_progress(pconf, 100) as meter1:
				meter1.moveto(100)
			with get_progress(pconf, 100) as meter2:
				meter2.moveto(100)
