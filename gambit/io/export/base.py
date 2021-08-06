from abc import ABC, abstractmethod
from typing import IO, Union

from gambit.query import QueryResults
from gambit.io import FilePath


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
