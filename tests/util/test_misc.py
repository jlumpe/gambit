import pytest
from string import ascii_letters

from gambit.util import misc


class TestZipStrict:
	"""Test the zip_strict() function."""
	N = 4

	def make_iterables(self, l):
		return [
			range(l),
			list(range(l)),
			iter(range(l)),  # IteraTOR, not IteraBLE
			ascii_letters[:l],
		]

	@pytest.mark.parametrize('l', [0, 10])
	@pytest.mark.parametrize('n', range(N))
	def test_valid(self, l, n):
		# Test with varying numbers of arguments
		iterables1 = self.make_iterables(l)[:n]
		iterables2 = self.make_iterables(l)[:n]
		assert list(misc.zip_strict(*iterables1)) == list(zip(*iterables2))

	@pytest.mark.parametrize('n', list(range(1, N)))
	@pytest.mark.parametrize('longer', [True, False])
	def test_invalid(self, n, longer):
		l = 10
		l2 = l + 1 if longer else l - 1

		for i in range(n+1):
			iterables = self.make_iterables(l)
			iterables.insert(i, range(l2))
			z = misc.zip_strict(*iterables)

			with pytest.raises(ValueError):
				list(z)


def test_chunk_slices():
	"""Test the chunk_slices() function."""

	for n, size in [(100, 10), (100, 30), (100, 1), (100, 1000)]:
		slices = list(misc.chunk_slices(n, size))
		ns = len(slices)

		for i, s in enumerate(slices):
			assert s.start == 0 if i == 0 else slices[i-1].stop
			assert s.step is None
			assert s.stop == s.start + size if i < ns - 1 else n

	assert list(misc.chunk_slices(0, 10)) == []


def test_join_list_human():
	l = ['foo', 'bar', 'baz']
	assert misc.join_list_human(l[:1]) == 'foo'
	assert misc.join_list_human(l[:2]) == 'foo and bar'
	assert misc.join_list_human(l[:3]) == 'foo, bar, and baz'
	assert misc.join_list_human(l[:1], 'or') == 'foo'
	assert misc.join_list_human(l[:2], 'or') == 'foo or bar'
	assert misc.join_list_human(l[:3], 'or') == 'foo, bar, or baz'
