"""Export results to JSON."""

import json
from typing import Union, IO, Any

from attr import attrs, attrib, asdict
from sqlalchemy.orm import Session

from gambit.query import QueryResultItem, QueryResults
from gambit.classify import ClassifierResult, GenomeMatch
from gambit.db import ReferenceGenomeSet, Taxon, AnnotatedGenome, Genome
import gambit.io.json as gjson
from gambit.io.util import FilePath, maybe_open
from gambit.util.misc import singledispatchmethod, type_singledispatchmethod
from gambit.util.typing import is_optional, unwrap_optional
from .base import asdict_default, BaseJSONResultsExporter


def _todict(obj, attrs):
	return {a: getattr(obj, a) for a in attrs}


@attrs()
class ResultsArchiveWriter(BaseJSONResultsExporter):
	"""Exports query results to "archive" format which captures all stored data.

	This format is not intended to be read by users of the application.
	The exported data can be read and converted back into an identical :class:`QueryResults`
	object using :class:`.ResultsArchiveReader`.

	Attributes
	----------
	install_info
		Add results of :func:`gambit.util.dev.install_info` to the ``QueryResults.extra`` dict.
	"""
	install_info: bool = attrib(default=False)

	to_json = singledispatchmethod(BaseJSONResultsExporter.to_json)

	to_json.register(ClassifierResult, asdict_default)
	to_json.register(GenomeMatch, asdict_default)
	to_json.register(QueryResultItem, asdict_default)

	@to_json.register(QueryResults)
	def _queryresults_to_json(self, results):
		data = asdict(results)

		if self.install_info:
			from gambit.util.dev import install_info
			data['extra']['install_info'] = install_info()

		return data

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

	@type_singledispatchmethod
	def _from_json(self, cls, data, ctx):
		"""Default implementation."""
		if is_optional(cls):
			if data is None:
				return None
			else:
				return self._from_json(unwrap_optional(cls), data, ctx)

		if hasattr(cls, '__attrs_attrs__'):
			return self._attrs_from_json(cls, data, ctx)
		else:
			return gjson.from_json(data, cls)

	def _attrs_from_json(self, cls, data, ctx, values=None):
		"""Create an attrs class instance from JSON data.

		``values`` is a dictionary of already-deserialized attribute values.
		"""
		kw = dict()

		for a in cls.__attrs_attrs__:
			if values is not None and a.name in values:
				kw[a.name] = values[a.name]
			else:
				atype = Any if a.type is None else a.type
				kw[a.name] = self._from_json(atype, data[a.name], ctx)

		return cls(**kw)

	@_from_json.register(ReferenceGenomeSet)
	def _genomeset_from_json(self, cls, data, ctx):
		assert data is not None
		return self.session.query(ReferenceGenomeSet).filter_by(key=data['key'], version=data['version']).one()

	@_from_json.register(AnnotatedGenome)
	def _genome_from_json(self, cls, data, ctx):
		key = data['key']
		gset_id = ctx['genomeset_id']
		return self.session.query(AnnotatedGenome)\
			.join(Genome)\
			.filter(AnnotatedGenome.genome_set_id == gset_id, Genome.key == key)\
			.one()

	@_from_json.register(Taxon)
	def _taxon_from_json(self, cls, data, ctx):
		key = data['key']
		gset_id = ctx['genomeset_id']
		return self.session.query(Taxon).filter_by(genome_set_id=gset_id, key=key).one()

	def results_from_json(self, data):
		genomeset = self._from_json(ReferenceGenomeSet, data['genomeset'], dict())

		# Add genome set to context so the correct AnnotatedGenomes can be loaded.
		ctx = dict(genomeset_id=genomeset.id)

		items = [self._from_json(QueryResultItem, item, ctx) for item in data['items']]
		return self._attrs_from_json(QueryResults, data, ctx, dict(genomeset=genomeset, items=items))

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
