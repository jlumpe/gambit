"""Export results to JSON."""

import json
import sys

from attr import attrs, attrib, asdict

from .base import AbstractResultsExporter
from gambit.query.results import QueryResults, QueryResultItem
from gambit.db.models import ReferenceGenomeSet, Taxon, AnnotatedGenome
import gambit.io.json as gjson


if sys.version_info[1] >= 8:
	from functools import singledispatchmethod

else:
	# Not available in 3.7, make simple implementation
	from functools import singledispatch, wraps

	def singledispatchmethod(func):
		dispatcher = singledispatch(func)

		@wraps(func)
		def wrapper(self, arg, *rest, **kw):
			impl = dispatcher.dispatch(type(arg))
			return impl(self, arg, *rest, **kw)

		wrapper.register = dispatcher.register
		return wrapper


def _todict(obj, attrs):
	return {a: getattr(obj, a) for a in attrs}


@attrs()
class JSONResultsExporter(AbstractResultsExporter):
	"""Exports query results in JSON format.

	Attributes
	----------
	dense
		Write with no whitespace to cut down on file size. Disable to produce more human-friendly
		output. Defaults to True.
	"""
	dense: bool = attrib(default=True)

	@singledispatchmethod
	def _to_json(self, obj):
		# Base case
		return gjson.to_json(obj)

	@_to_json.register(QueryResults)
	def _results_to_json(self, results: QueryResults):
		return asdict(results)

	@_to_json.register(QueryResultItem)
	def _item_to_json(self, item: QueryResultItem):
		return asdict(item)

	@_to_json.register(ReferenceGenomeSet)
	def _genomeset_to_json(self, gset: ReferenceGenomeSet):
		return _todict(gset, ['id', 'key', 'version', 'name', 'description'])

	@_to_json.register(Taxon)
	def _taxon_to_json(self, taxon: Taxon):
		return _todict(taxon, ['id', 'key', 'name', 'ncbi_id', 'distance_threshold'])

	@_to_json.register(AnnotatedGenome)
	def _genome_to_json(self, genome: AnnotatedGenome):
		data = _todict(genome, ['key', 'description', 'organism', 'taxon', 'ncbi_db', 'ncbi_id', 'genbank_acc', 'refseq_acc'])
		data['id'] = genome.genome_id
		return data

	def export(self, f, results: QueryResults):
		json.dump(results, f, default=self._to_json)
