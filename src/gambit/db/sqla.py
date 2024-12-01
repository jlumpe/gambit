"""Custom types and other utilities for SQLAlchemy."""
import os
from typing import Optional

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


def default_sessionmaker(bind, *, readonly: bool = True, class_: Optional[type] = None, **kw) -> sessionmaker:
	"""Create an SQLAlchemy ``sessionmaker`` using some common default settings.

	Parameters
	----------
	bind
		First argument to :class:`sqlalchemy.orm.sessionmaker`.
	readonly
		Sets the default value for the ``class_`` keyword argument (:class:`.ReadOnlySession` if True,
		otherwise uses the standard SQLAlchemy session type).
	\\**kw
		Additional keyword arguments to :class:`sqlalchemy.orm.sessionmaker`.
	"""
	if class_ is None:
		class_ = ReadOnlySession if readonly else Session
	# future=True - forwards compatibility with SQLAlchemy 2.0
	return sessionmaker(bind, class_=class_, future=True, **kw)


def file_sessionmaker(path: 'FilePath', **kw) -> sessionmaker:
	"""Get an SQLAlchemy ``sessionmaker`` for an sqlite database file.

	Parameters
	----------
	path
		Path to database file.
	\\**kw
		Additional keyword arguments to :func:`.default_sessionmaker` /
		:class:`sqlalchemy.orm.sessionmaker`.
	"""
	engine = create_engine(f'sqlite:///{os.fspath(path)}')
	return default_sessionmaker(engine, **kw)
