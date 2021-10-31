"""Export query results to CSV."""

import csv
from typing import Dict, Any, List, Union, Iterable, TextIO

from .base import AbstractResultsExporter
from gambit.query import QueryResultItem, QueryResults
from gambit.io.util import FilePath, maybe_open


def getattr_nested(obj, attrs: Union[str, Iterable[str]], pass_none=False):
	if isinstance(attrs, str):
		attrs = attrs.split('.')

	for attr in attrs:
		if pass_none and obj is None:
			return None

		obj = getattr(obj, attr)

	return obj


class CSVResultsExporter(AbstractResultsExporter):
	"""Exports query results in CSV format.

	Attributes
	----------
	format_opts
		Dialect and other formatting arguments passed to :func:`csv.write`.
	"""
	format_opts: Dict[str, Any]

	COLUMNS = [
		('query.name', 'input.label'),
		('query.path', 'input.file.path'),
		('predicted.name', 'report_taxon.name'),
		('predicted.rank', 'report_taxon.rank'),
		('predicted.ncbi_id', 'report_taxon.ncbi_id'),
		('predicted.threshold', 'report_taxon.distance_threshold'),
		('closest.distance', 'classifier_result.closest_match.distance'),
		('closest.description', 'classifier_result.closest_match.genome.description'),
	]

	def __init__(self, **format_opts):
		if 'dialect' not in format_opts:
			format_opts.setdefault('lineterminator', '\n')
			format_opts.setdefault('quoting', csv.QUOTE_MINIMAL)
		self.format_opts = format_opts

	def get_header(self) -> List[str]:
		"""Get values for header row."""
		return [name for name, _ in self.COLUMNS]

	def get_row(self, item: QueryResultItem) -> List:
		"""Get row values for single result item."""
		return [getattr_nested(item, attrs, pass_none=True) for _, attrs in self.COLUMNS]

	def export(self, file_or_path: Union[FilePath, TextIO], results: QueryResults):
		with maybe_open(file_or_path, 'w') as f:
			writer = csv.writer(f, **self.format_opts)

			writer.writerow(self.get_header())
			for item in results.items:
				writer.writerow(self.get_row(item))
