"""Export query results in various formats."""

import json
from abc import ABC, abstractmethod
from typing import IO, Union, TextIO, Any, Iterable
import csv
from functools import singledispatchmethod

from attr import attrs, attrib, asdict
from sqlalchemy.orm import Session

from gambit.util.io import FilePath, maybe_open
import gambit.util.json as gjson
from gambit.query import QueryResults, QueryResultItem, QueryInput
from gambit.db import ReferenceGenomeSet, Taxon, AnnotatedGenome, Genome


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
		Dialect and other formatting arguments passed to :func:`csv.writer`.
	"""
	format_opts: dict[str, Any]

	COLUMNS = [
		('query', 'input.label'),
		('predicted.name', 'report_taxon.name'),
		('predicted.rank', 'report_taxon.rank'),
		('predicted.ncbi_id', 'report_taxon.ncbi_id'),
		('predicted.threshold', 'report_taxon.distance_threshold'),
		('closest.distance', 'classifier_result.closest_match.distance'),
		('closest.description', 'classifier_result.closest_match.genome.description'),
		('next.name', 'classifier_result.next_taxon.name'),
		('next.rank', 'classifier_result.next_taxon.rank'),
		('next.ncbi_id', 'classifier_result.next_taxon.ncbi_id'),
		('next.threshold', 'classifier_result.next_taxon.distance_threshold'),
	]

	def __init__(self, **format_opts):
		if 'dialect' not in format_opts:
			format_opts.setdefault('lineterminator', '\n')
			format_opts.setdefault('quoting', csv.QUOTE_MINIMAL)
		self.format_opts = format_opts

	def get_header(self) -> list[str]:
		"""Get values for header row."""
		return [name for name, _ in self.COLUMNS]

	def get_row(self, item: QueryResultItem) -> list:
		"""Get row values for single result item."""
		return [getattr_nested(item, attrs, pass_none=True) for _, attrs in self.COLUMNS]

	def export(self, file_or_path: Union[FilePath, TextIO], results: QueryResults):
		with maybe_open(file_or_path, 'w') as f:
			writer = csv.writer(f, **self.format_opts)

			writer.writerow(self.get_header())
			for item in results.items:
				writer.writerow(self.get_row(item))


@attrs()
class JSONResultsExporter(BaseJSONResultsExporter):
	"""Exports query results in basic JSON format.

	Currently it assumes that the query was run with ``classify_strict=False``, so the only
	relevant information from ``ClassifierResult`` is the closest genome match.
	"""

	to_json = singledispatchmethod(BaseJSONResultsExporter.to_json)

	@to_json.register(QueryResults)
	def _results_to_json(self, results: QueryResults):
		data = asdict(results, recurse=False)
		del data['params']  # Parameters not currently exposed thru CLI, so omit for now.
		return data

	@to_json.register(QueryResultItem)
	def _item_to_json(self, item: QueryResultItem):
		return dict(
			query=item.input,
			predicted_taxon=item.report_taxon,
			next_taxon=item.classifier_result.next_taxon,
			closest_genomes=item.closest_genomes,
		)

	@to_json.register(QueryInput)
	def _input_to_json(self, input: QueryInput):
		return dict(
			name=input.label,
			path=None if input.file is None else input.file.path,
			format=None if input.file is None else input.file.format,
		)

	@to_json.register(ReferenceGenomeSet)
	def _genomeset_to_json(self, gset: ReferenceGenomeSet):
		return _todict(gset, ['id', 'key', 'version', 'name', 'description'])

	@to_json.register(Taxon)
	def _taxon_to_json(self, taxon: Taxon):
		return _todict(taxon, ['id', 'key', 'name', 'ncbi_id', 'rank', 'distance_threshold'])

	@to_json.register(AnnotatedGenome)
	def _genome_to_json(self, genome: AnnotatedGenome):
		data = _todict(genome, ['key', 'description', 'organism', 'ncbi_db', 'ncbi_id', 'genbank_acc', 'refseq_acc'])
		data['id'] = genome.genome_id
		data['taxonomy'] = list(genome.taxon.ancestors(incself=True))
		return data


class ResultsArchiveWriter(BaseJSONResultsExporter):
	"""Exports query results to "archive" format which captures all stored data.

	This format is not intended to be read by users of the application.
	The exported data can be read and converted back into an identical
	:class:`~gambit.query.QueryResults` object using :class:`.ResultsArchiveReader`.

	Only the ID attributes of database models are saved, when loading the saved results the models
	are recreated by database queries.
	"""

	to_json = singledispatchmethod(BaseJSONResultsExporter.to_json)

	@to_json.register(ReferenceGenomeSet)
	def _genomeset_to_json(self, gset: ReferenceGenomeSet):
		return _todict(gset, ['key', 'version'])

	@to_json.register(Taxon)
	def _taxon_to_json(self, taxon: Taxon):
		return _todict(taxon, ['key'])

	@to_json.register(AnnotatedGenome)
	def _genome_to_json(self, genome: AnnotatedGenome):
		return _todict(genome, ['key'])


class ResultsArchiveReader:
	"""Loads query results from file created by :class:`ResultsArchiveWriter`.

	Attributes
	----------
	session
		SQLAlchemy session used to load database objects.
	"""
	session: Session

	def __init__(self, session):
		self.session = session

		self._init_converter()

		# Loading the Taxon and AnnotatedGenome instances from the database requires not just their
		# ID (key attribute) values but also the ReferenceGenomeSet they belong to. Setting this
		# attribute to the genome set instance of the results currently being loaded is a somewhat
		# hacky method of passing this information to the unstructuring hook functions. There isn't
		# a much better way of doing this without reimplementing a lot of the cattrs machinery.
		self._current_genomeset = None

	def _init_converter(self):
		"""Initialize the cattrs converter instance.

		This is a clone of the converter instance in gambit.util.json, with additional structuring
		hooks registered to methods on this instance.
		"""
		self._converter = gjson.converter.copy()
		self._converter.register_structure_hook(ReferenceGenomeSet, self._structure_genomeset)
		self._converter.register_structure_hook(AnnotatedGenome, self._structure_genome)
		self._converter.register_structure_hook(Taxon, self._structure_taxon)

	def read(self, file_or_path: Union[FilePath, IO]) -> QueryResults:
		"""Read query results from JSON file.

		Parameters
		----------
		file_or_path
			Readable file object or file path.
		"""
		with maybe_open(file_or_path) as f:
			data = json.load(f)

		return self.results_from_json(data)

	def results_from_json(self, data: dict[str, Any]) -> QueryResults:
		"""Recreate results object from loaded JSON data."""

		gset_key = data['genomeset']['key']
		gset_version = data['genomeset']['version']
		self._current_genomeset =  self.session.query(ReferenceGenomeSet) \
			.filter_by(key=gset_key, version=gset_version) \
			.one()

		try:
			return self._converter.structure(data, QueryResults)

		finally:
			self._current_genomeset = None

	def _structure_genomeset(self, data: dict[str, Any], cls=None):
		return self._current_genomeset

	def _structure_genome(self, data: dict[str, Any], cls=None) -> AnnotatedGenome:
		key = data['key']
		gset_id = self._current_genomeset.id
		return self.session.query(AnnotatedGenome)\
			.join(Genome)\
			.filter(AnnotatedGenome.genome_set_id == gset_id, Genome.key == key)\
			.one()

	def _structure_taxon(self, data: dict[str, Any], cls=None) -> Taxon:
		key = data['key']
		gset_id = self._current_genomeset.id
		return self.session.query(Taxon).filter_by(genome_set_id=gset_id, key=key).one()
