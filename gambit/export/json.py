"""Export results to JSON."""

import json
from typing import Union, IO

from attr import attrs, attrib, asdict

from .base import AbstractResultsExporter
from gambit.query import QueryResultItem, QueryResults
from gambit.classify import ClassifierResult
from gambit.db import ReferenceGenomeSet, Taxon, AnnotatedGenome
import gambit.io.json as gjson
from gambit.io.util import FilePath, maybe_open
from gambit.util.misc import singledispatchmethod


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

	def export(self, file_or_path: Union[FilePath, IO], results: QueryResults):
		opts = dict(indent=4, sort_keys=True) if self.pretty else dict()
		with maybe_open(file_or_path, 'w') as f:
			json.dump(results, f, default=self.to_json, **opts)


@attrs()
class JSONResultsExporter(BaseJSONResultsExporter):
	"""Exports query results in JSON format.
	"""

	to_json = singledispatchmethod(BaseJSONResultsExporter.to_json)

	to_json.register(ClassifierResult, asdict_default)
	to_json.register(QueryResults, asdict_default)
	to_json.register(QueryResultItem, asdict_default)

	@to_json.register(ReferenceGenomeSet)
	def _genomeset_to_json(self, gset: ReferenceGenomeSet):
		return _todict(gset, ['id', 'key', 'version', 'name', 'description'])

	@to_json.register(Taxon)
	def _taxon_to_json(self, taxon: Taxon):
		return _todict(taxon, ['id', 'key', 'name', 'ncbi_id', 'distance_threshold'])

	@to_json.register(AnnotatedGenome)
	def _genome_to_json(self, genome: AnnotatedGenome):
		data = _todict(genome, ['key', 'description', 'organism', 'taxon', 'ncbi_db', 'ncbi_id', 'genbank_acc', 'refseq_acc'])
		data['id'] = genome.genome_id
		return data
