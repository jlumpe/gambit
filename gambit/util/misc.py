"""Utility code that doesn't fit anywhere else."""

import sys
from typing import Iterator, Tuple


def zip_strict(*iterables: Iterator) -> Iterator[Tuple]:
	"""Like the builtin ``zip`` function but raises an error if any argument is exhausted before the others.

	Parameters
	----------
	iterables
		Any number of iterable objects.

	Raises
	------
	ValueError
	"""
	# Builtin zip gives empty output on empty input
	if not iterables:
		return

	itrs = list(map(iter, iterables))
	n = len(itrs)

	while True:
		# Try to get next item of first iterator
		try:
			first = next(itrs[0])
		except StopIteration:
			# First ran out, make sure others do too
			for i in range(1, n):
				try:
					next(itrs[i])
				except StopIteration:
					pass
				else:
					raise ValueError(f'Iterable {i} yielded more items than previous iterables')

			return

		# First iterator not done, get rest
		out = [first]
		for i in range(1, n):
			try:
				value = next(itrs[i])
			except StopIteration:
				raise ValueError(f'Iterable {i} exhausted before previous iterables')

			out.append(value)

		yield tuple(out)


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


if sys.version_info[1] >= 8:
	from functools import singledispatchmethod

else:
	# Not available in 3.7, make simple implementation
	from functools import singledispatch, wraps

	def singledispatchmethod(func):
		dispatcher = singledispatch(func)

		@wraps(func)
		def wrapper(self, arg, *rest, **kw):
			impl = dispatcher.dispatch(type(arg))
			return impl(self, arg, *rest, **kw)

		wrapper.register = dispatcher.register
		return wrapper
