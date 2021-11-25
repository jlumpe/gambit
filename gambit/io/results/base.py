import json
from abc import ABC, abstractmethod
from typing import IO, Union, TextIO

from attr import asdict, attrs, attrib

from gambit.io.util import FilePath, maybe_open
import gambit.io.json as gjson
from gambit.query import QueryResults


class AbstractResultsExporter(ABC):
	"""Base for classes that export formatted query results.

	Subclasses must implement :meth:`export`.
	"""

	@abstractmethod
	def export(self, file_or_path: Union[FilePath, IO], results: QueryResults):
		"""Write query results to file.

		Parameters
		----------
		file_or_path
			Open file-like object or file path to write to.
		results
			Results to export.
		"""


def _todict(obj, attrs):
	return {a: getattr(obj, a) for a in attrs}


def asdict_method(recurse=False, **kw):
	"""Create a ``to_json`` method which calls :func:`attrs.asdict` with the given options."""
	def method(self, obj):
		return asdict(obj, recurse=recurse, **kw)
	return method


asdict_default = asdict_method()


@attrs()
class BaseJSONResultsExporter(AbstractResultsExporter):
	"""Base class for JSON exporters.

	Subclasses need to implement the ``to_json`` method.

	Attributes
	----------
	pretty
		Write in more human-readable but less compact format. Defaults to False.
	"""
	pretty: bool = attrib(default=False)

	def to_json(self, obj):
		"""Convert object to JSON-compatible format (need not work recursively)."""
		return gjson.to_json(obj)

	def export(self, file_or_path: Union[FilePath, TextIO], results: QueryResults):
		opts = dict(indent=4, sort_keys=True) if self.pretty else dict()
		with maybe_open(file_or_path, 'w') as f:
			json.dump(results, f, default=self.to_json, **opts)
