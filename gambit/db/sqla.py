"""Custom types and other utilities for SQLAlchemy."""

import json

from sqlalchemy.types import TypeDecorator, String
from sqlalchemy.orm.session import Session

import gambit.io.json as gjson


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
