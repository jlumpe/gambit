"""Utility code that doesn't fit anywhere else."""

import sys
from typing import Iterator, Tuple, Callable, Iterable
from functools import singledispatch, wraps


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


# singledispatchmethod ot available in 3.7
if sys.version_info[1] >= 8:
	from functools import singledispatchmethod

else:
	# Make simple implementation
	def singledispatchmethod(func):
		dispatcher = singledispatch(func)

		@wraps(func)
		def wrapper(self, arg, *rest, **kw):
			impl = dispatcher.dispatch(type(arg))
			return impl(self, arg, *rest, **kw)

		wrapper.register = dispatcher.register
		return wrapper


def type_singledispatchmethod(func: Callable):
	"""
	Similar to ``singledispatchmethod``, but the first (non-self) argument is expected to be a
	type and dispatch occurs on the argument's *value*.

	Parameters
	----------
	func
		Default implementation. Signature should start with ``(self, cls: type, ...)``.

	Returns
	-------
	Callable
		Function with ``register`` and ``dispatch`` attributes similar to ``singledispatchmethod``.
	"""

	dispatcher = singledispatch(func)

	@wraps(func)
	def wrapper(self, cls, *rest, **kw):
		if isinstance(cls, type):
			impl = dispatcher.dispatch(cls)
		else:
			# Use default implementation for non-types (e.g. stuff from tye typing module)
			impl = dispatcher.dispatch(object)

		return impl(self, cls, *rest, **kw)

	wrapper.register = dispatcher.register
	wrapper.dispatch = dispatcher.dispatch
	return wrapper


def is_importable(module: str) -> bool:
	"""Check if the specified module is importable, without actually importing it."""
	from importlib.util import find_spec
	return find_spec(module) is not None


def join_list_human(strings: Iterable[str], conj: str='and') -> str:
	"""Join items into a single human-readable string with commas and the given conjunction."""
	strings = list(strings)
	if len(strings) > 2:
		s = ', '.join(strings[:-1]) + ','
		return ' '.join([s, conj, strings[-1]])
	if len(strings) == 2:
		return ' '.join([strings[0], conj, strings[1]])
	if len(strings) == 1:
		return strings[0]
	return ''
