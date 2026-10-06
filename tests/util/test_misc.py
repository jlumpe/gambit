from gambit.util import misc


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
