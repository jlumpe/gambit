from gambit.util.misc import chunk_slices


def test_chunk_slices():
	"""Test the chunk_slices() function."""

	for n, size in [(100, 10), (100, 30), (100, 1), (100, 1000)]:
		slices = list(chunk_slices(n, size))
		ns = len(slices)

		for i, s in enumerate(slices):
			assert s.start == 0 if i == 0 else slices[i-1].stop
			assert s.step is None
			assert s.stop == s.start + size if i < ns - 1 else n

	assert list(chunk_slices(0, 10)) == []
