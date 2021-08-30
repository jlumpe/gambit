"""Export results to JSON."""

import json
from typing import Union, TextIO

from attr import attrs, attrib, asdict

from .base import AbstractResultsExporter
from gambit.query import QueryResultItem, QueryResults, QueryInput
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

	def export(self, file_or_path: Union[FilePath, TextIO], results: QueryResults):
		opts = dict(indent=4, sort_keys=True) if self.pretty else dict()
		with maybe_open(file_or_path, 'w') as f:
			json.dump(results, f, default=self.to_json, **opts)


@attrs()
class JSONResultsExporter(BaseJSONResultsExporter):
	"""Exports query results in basic JSON format.

	Currently it assumes that the query was run with ``classify_strict=False``\\ , so the only
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
			closest_genome=item.classifier_result.closest_match.genome,
			closest_genome_distance=item.classifier_result.closest_match.distance,
		)

	@to_json.register(QueryInput)
	def _input_to_json(self, input: QueryInput):
		return dict(name=input.label, path=input.file.path, format=input.file.format)

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
