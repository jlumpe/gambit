"""Utility code for reading/writing data files."""

import os
from typing import Union, Optional

#: Alias for types which can represent a file system path
FilePath = Union[str, os.PathLike]


COMPRESSED_OPENERS = {None: open}


def _compressed_opener(compression):
	"""Decorator to register opener functions for compression types."""
	def decorator(func):
		COMPRESSED_OPENERS[compression] = func
		return func
	return decorator


@_compressed_opener('gzip')
def _open_gzip(path, mode, **kwargs):
	"""Opener for gzip-compressed files."""
	import gzip

	if mode is None:
		mode = 'rt'

	# gzip defaults to binary mode, change to text instead of not specified
	if mode[-1] not in 'tb':
		mode += 't'

	return gzip.open(path, mode=mode, **kwargs)


def open_compressed(compression: Optional[str], path: FilePath, mode: Optional[str] = None, **kwargs):
	"""Open a file with compression method specified by a string.

	Parameters
	----------
	compression : str
		Compression method. None is no compression. Keys of :data:`COMPRESSED_OPENERS` are the
		allowed values.
	path
		Path of file to open. May be string or path-like object.
	mode : str
		Mode to open file in - same as in :func:`open`.
	\\**kwargs
		Additional text-specific keyword arguments identical to the following :func:`open`
		arguments: ``encoding``, ``errors``, and ``newlines``.

	Returns
	-------
	io.BufferedIOBase
		Open file object.
	"""

	try:
		opener = COMPRESSED_OPENERS[compression]

	except KeyError:
		raise ValueError(f'Unknown compression type {compression!r}') from None

	return opener(os.fsdecode(path), mode=mode, **kwargs)


class ClosingIterator:
	"""Wraps an iterator which reads from a stream, closes the stream when finished.

	Used to wrap return values from functions which do some sort of lazy IO
	operation (specifically :func:`Bio.SeqIO.parse`) and return an iterator
	which reads from a stream every time ``next()`` is called on it. The object
	is an iterator itself, but will close the stream automatically when it
	finishes. May also be used as a context manager which closes the stream
	on exit.

	Attributes
	----------
	fobj

		The underlying file-like object or stream which the instance is
		responsible for closing

	iterator

		The iterator which the instance wraps.

	closed

		Read-only boolean property, mirrors the same attribute of :attr:`fobj`.

	Parameters
	----------
	iterable
		Iterable to iterate over. The :attr:`iterator` attribute will be obtained from calling
		:func:`iter` on this.
	fobj
		File-like object to close when the iterator finishes, context is exited or the :meth:`close`
		method is called.
	"""

	def __init__(self, iterable, fobj):
		self.iterator = iter(iterable)
		self.fobj = fobj

	def __iter__(self):
		return self

	def __next__(self):
		try:
			return next(self.iterator)

		except StopIteration:
			# Close when iterator runs out
			self.close()
			raise

	def close(self):
		"""Close the stream.

		Just calls the ``close`` method on :attr:`fobj`.
		"""
		self.fobj.close()

	@property
	def closed(self):
		return self.fobj.closed

	def __enter__(self):
		return self

	def __exit__(self, *args):
		self.close()
