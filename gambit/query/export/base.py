from abc import ABC, abstractmethod
from gambit.query.results import QueryResults


class AbstractResultsExporter(ABC):
	"""Base for classes that export formatted query results.

	Subclasses must implement :meth:`export`.
	"""

	@abstractmethod
	def export(self, f, results: QueryResults):
		"""Write query results to file.

		Parameters
		----------
		f
			Writable file-like object.
		results
			Results to export.
		"""
