"""Export results to JSON."""

import json
from typing import Union, IO, Any
from functools import singledispatchmethod

from sqlalchemy.orm import Session

from gambit.query import QueryResults
from gambit.db import ReferenceGenomeSet, Taxon, AnnotatedGenome, Genome
import gambit.util.json as gjson
from gambit.util.io import FilePath, maybe_open
from .base import BaseJSONResultsExporter, _todict


class ResultsArchiveWriter(BaseJSONResultsExporter):
	"""Exports query results to "archive" format which captures all stored data.

	This format is not intended to be read by users of the application.
	The exported data can be read and converted back into an identical :class:`QueryResults`
	object using :class:`.ResultsArchiveReader`.

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
