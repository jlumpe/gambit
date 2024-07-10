"""Utility code for reading/writing data files."""

import os
from io import TextIOWrapper
from typing import Union, Optional, IO, BinaryIO, ContextManager, Iterable, TypeVar
from contextlib import nullcontext

#: Alias for types which can represent a file system path
FilePath = Union[str, os.PathLike]

T = TypeVar('T')

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
	return gzip.open(path, mode=mode, **kwargs)


@_compressed_opener('auto')
def _open_auto(path, mode, **kwargs):
	"""Open file for reading with compression determined automatically."""

	if mode[0] != 'r':
		raise ValueError('Automatic compression detection only supported for reading.')

	file = open(path, 'rb')

	try:
		compression = guess_compression(file)
		file.seek(0)

		if compression is None:
			binary = file
		elif compression == 'gzip':
			import gzip
			binary = gzip.GzipFile(fileobj=file, mode='rb')
		else:
			assert 0

		return TextIOWrapper(binary, **kwargs) if mode[1] == 't' else binary

	except Exception:
		file.close()


def guess_compression(fobj: BinaryIO) -> Optional[str]:
	"""Guess the compression mode of an readable file-like object in binary mode.

	Assumes the current position is at the beginning of the file.
	"""
	magic = fobj.read(2)

	if magic == b'\x1f\x8b':
		return 'gzip'
	else:
		return None


def open_compressed(compression: Optional[str],
                    path: FilePath,
                    mode: str = 'rt',
                    **kwargs,
                    ) -> IO:
	"""Open a file with compression method specified by a string.

	Parameters
	----------
	compression : str
		Compression method. None is no compression. Keys of :data:`COMPRESSED_OPENERS` are the
		allowed values.
	path
		Path of file to open. May be string or path-like object.
	mode : str
		Mode to open file in - similar to :func:`open`. Must be exactly two characters, the first
		in ``rwax`` and the second in``tb``.
	\\**kwargs
		Additional text-specific keyword arguments identical to the following :func:`open`
		arguments: ``encoding``, ``errors``, and ``newlines``.

	Returns
	-------
	IO
		Open file object.
	"""
	if not(len(mode) == 2 and mode[0] in 'rwax' and mode[1] in 'tb'):
		msg = f'Invalid mode {mode!r}'
		if mode in 'rwax':
			msg += ' (must specify either binary or text mode)'
		raise ValueError(msg)

	try:
		opener = COMPRESSED_OPENERS[compression]

	except KeyError:
		raise ValueError(f'Unknown compression type {compression!r}') from None

	return opener(os.fsdecode(path), mode=mode, **kwargs)


class ClosingIterator(Iterable[T]):
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
	def closed(self) -> bool:
		return self.fobj.closed

	def __enter__(self):
		return self

	def __exit__(self, *args):
		self.close()


def maybe_open(file_or_path: Union[FilePath, IO], mode: str = 'r', **open_kw) -> ContextManager[IO]:
	"""Open a file given a file path as an argument, but pass existing file objects though.

	Intended to be used by API functions that take either type as an argument. If a file path is
	passed the function will need to call ``open`` to get the file object to use, and will need to
	close that object after it is done. If an existing file object is passed, it should be left to
	the caller of the function to close it afterwards. This function returns a context manager which
	performs the correct action for both opening and closing.

	Parameters
	----------
	file_or_path
		A path-like object or open file object.
	mode
		Mode to open file in.
	\\**open_kw
		Keyword arguments to :func:`open`.

	Returns
	-------
	ContextManager[IO]
		Context manager which gives an open file object on enter and closes it on exit only if it
		was opened by this function.
	"""
	try:
		# Try to interpret as path
		path = os.fspath(file_or_path)
	except TypeError:
		# Not a path, assume file object
		# Return context manager which gives this object on enter and does not close on exit
		return nullcontext(file_or_path)
	else:
		# Is a path, just open
		return open(path, mode, **open_kw)


def read_lines(file_or_path: Union[FilePath, IO], strip: bool=True, skip_empty: bool=False) -> Iterable[str]:
	"""Iterate over lines in text file.

	Parameters
	----------
	file_or_path
		A path-like object or open file object.
	strip
		Strip whitespace from lines.
	skip_empty
		Omit empty lines.

	Returns
	-------
	Iterable[str]
		Iterator over lines with trailing newlines removed.
	"""
	with maybe_open(file_or_path) as file:
		for line in file:
			line = line.strip() if strip else line.rstrip('\n')
			if not (skip_empty and not line):
				yield line


def write_lines(lines: Iterable, file_or_path: Union[FilePath, IO]):
	"""Write strings to text file, one per line.

	Parameters
	----------
	lines
		Iterable of lines to write.
	file_or_path
		A path-like object or open file object.
	"""
	with maybe_open(file_or_path, 'w') as file:
		file.writelines(str(line) + '\n' for line in lines)
