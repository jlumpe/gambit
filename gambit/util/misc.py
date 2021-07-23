"""Utility code that doesn't fit anywhere else."""

from typing import Iterator


def chunk_slices(n: int, size: int) -> Iterator[slice]:
	"""Iterate over slice objects which split a sequence of length ``n`` into chunks of size ``size``.

	Parameters
	----------
	n
		Length of sequence.
	size
		Size of chunks (apart from last).
	"""
	if size <= 0:
		raise ValueError('Size must be positive')

	start = 0
	while start < n:
		stop = start + size
		yield slice(start, stop)
		start = stop
