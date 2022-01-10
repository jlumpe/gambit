"""Custom types and other utilities for SQLAlchemy."""
import os

from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker, Session
from sqlalchemy.types import TypeDecorator, String

import gambit.util.json as gjson
from gambit.util.io import FilePath


class ReadOnlySession(Session):
	"""Session class that doesn't allow flushing/committing."""

	def flush(self, *args, **kwargs):
		# Make flush a no-op
		pass

	def commit(self):
		raise TypeError('Session is read-only')


class JsonString(TypeDecorator):
	"""SQLA column type for JSON data which is stored in the database as a standard string column.

	Data is automatically serialized/unserialized when saved/loaded.
	Important: mutation tracking is not enabled for this type. If the value is a list or dict and
	you modify it in place these changes will not be detected. Instead, re-assign the attribute.
	"""

	impl = String

	def process_bind_param(self, value, dialect):
		return None if value is None else gjson.dumps(value)

	def process_result_value(self, value, dialect):
		return None if value is None else gjson.loads(value)


def file_sessionmaker(path: FilePath, readonly: bool = True, cls: type = None, **kw) -> sessionmaker:
	"""Get an SQLAlchemy ``sessionmaker`` for an sqlite database file.

	Parameters
	----------
	path
		Path to database file.
	readonly
		Sets the default value for ``class_``.
	cls
		SQLAlchemy ``Session`` subclass to use. Defaults to :class:`gambit.db.sqla.ReadOnlySession`
		if ``readonly=True``, otherwise uses the standard SQLAlchemy session type.
	\\**kw
		Additional keyword arguments to :class:`sqlalchemy.orm.sessionmaker`.
	"""
	if cls is None:
		cls = ReadOnlySession if readonly else Session
	engine = create_engine(f'sqlite:///{os.fspath(path)}')
	return sessionmaker(engine, class_=cls, **kw)
