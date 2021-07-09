"""Query a GAMBIT database."""

from typing import Sequence, Optional, Union

import numpy as np

from gambit.db.gambitdb import GAMBITDatabase
from gambit.kmers import KmerSignature
from gambit.io.seq import SequenceFile
from gambit.metric import jaccard_sparse_array
from .classify import find_matches, consensus_taxon, reportable_taxon, matching_taxon
from .results import QueryInput, GenomeMatch, QueryResultItem, QueryResults


def _taxon_repr(taxon):
	"""Get a short string representation of a Taxon for logging and warning/error messages."""
	return f'{taxon.id}:{taxon.name}'


def runquery(db: GAMBITDatabase,
             queries: Sequence[KmerSignature],
             inputs: Optional[Sequence[Union[QueryInput, SequenceFile, str]]],
             ) -> QueryResults:
	"""Predict the taxonomy of one or more query genomes using a given GAMBIT reference database.

	Parameters
	----------
	db
		Database to query.
	queries
		Sequence of k-mer signatures.
	inputs
		Description for each input, converted to :class:`gambit.query.result.QueryInput` in results
		object. Only used for reporting, does not any other aspect of results. Items can be
		``QueryInput``, ``SequenceFile`` or ``str``.
	"""
	queries = list(queries)

	if len(queries) == 0:
		raise ValueError('Must supply at least one query.')

	if inputs is not None:
		inputs = list(map(QueryInput.convert, inputs))
		if len(inputs) != len(queries):
			raise ValueError('Number of inputs does not match number of queries.')
	else:
		inputs = [QueryInput(str(i + 1)) for i in range(len(queries))]

	items = [_query_single(db, query, input) for query, input in zip(queries, inputs)]
	return QueryResults(items=items, genomeset=db.genomeset, signaturesmeta=db.signatures_meta)


def _query_single(db: GAMBITDatabase, sig: np.ndarray, input: QueryInput):
	dists = jaccard_sparse_array(sig, db.genome_signatures, distance=True)
	matches = find_matches(zip(db.genomes, dists))
	consensus, others = consensus_taxon(matches.keys())

	# Find closest match
	closest = np.argmin(dists)
	closest_match = GenomeMatch(
		genome=db.genomes[closest],
		distance=dists[closest],
		matching_taxon=matching_taxon(db.genomes[closest].taxon, dists[closest]),
	)

	# No matches found
	if not matches:
		return QueryResultItem(
			input=input,
			success=True,
			primary_match=None,
			closest_match=closest_match,
			predicted_taxon=None,
			report_taxon=None,
		)

	# Find primary match
	if consensus is None:
		primary_match = None

	else:
		best_i = None
		best_d = float('inf')
		best_taxon = None

		for taxon, idxs in matches.items():
			if consensus not in taxon.ancestors(incself=True):
				continue

			for i in idxs:
				if dists[i] < best_d:
					best_i = i
					best_d = dists[i]
					best_taxon = taxon

		assert best_i is not None
		primary_match = GenomeMatch(
			genome=db.genomes[best_i],
			distance=best_d,
			matching_taxon=best_taxon,
		)

	item = QueryResultItem(
		input=input,
		success=True,
		primary_match=primary_match,
		closest_match=closest_match,
		predicted_taxon=consensus,
		report_taxon=None if consensus is None else reportable_taxon(consensus),
	)

	# Warn of inconsistent matches
	if others:
		msg = f'Query matched {len(others)} inconsistent taxa: '
		msg += ', '.join(map(_taxon_repr, others))
		msg += '. Reporting lowest common ancestor of this set.'
		item.warnings.append(msg)

	# No consensus found - matches do not have common ancestor
	if consensus is None:
		item.success = False
		item.error = 'Matched taxa have no common ancestor.'

	# Consensus has no reportable ancestor
	elif item.report_taxon is None:
		item.success = False
		item.error = (
			f'Matched taxon {_taxon_repr(consensus)} has no reportable ancestor. '
			'This indicates a problem with the database.'
		)

	return item
