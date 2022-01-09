"""Export results to JSON."""

from attr import attrs, asdict

from .base import _todict, BaseJSONResultsExporter
from gambit.query import QueryResultItem, QueryResults, QueryInput
from gambit.db import ReferenceGenomeSet, Taxon, AnnotatedGenome
from gambit.util.misc import singledispatchmethod


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
			closest_genome=item.classifier_result.closest_match.genome,
			closest_genome_distance=item.classifier_result.closest_match.distance,
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
