"""Query a GAMBIT database."""

from typing import Sequence, Optional, Union

import numpy as np

from gambit.db.gambitdb import GAMBITDatabase
from gambit.kmers import KmerSignature
from gambit.io.seq import SequenceFile
from gambit.util.misc import zip_strict
from gambit.util.progress import progress_config, iter_progress
from gambit.metric import jaccard_sparse_matrix
from .classify import find_matches, consensus_taxon, reportable_taxon, matching_taxon
from .results import QueryInput, GenomeMatch, QueryResultItem, QueryResults


def _taxon_repr(taxon):
	"""Get a short string representation of a Taxon for logging and warning/error messages."""
	return f'{taxon.id}:{taxon.name}'


def runquery(db: GAMBITDatabase,
             queries: Sequence[KmerSignature],
             inputs: Optional[Sequence[Union[QueryInput, SequenceFile, str]]] = None,
             *,
             chunksize: Optional[int] = 1000,
             progress = None,
             ) -> QueryResults:
	"""Predict the taxonomy of one or more query genomes using a GAMBIT reference database.

	Parameters
	----------
	db
		Database to query.
	queries
		Sequence of k-mer signatures for query genomes.
	inputs
		Description for each input, converted to :class:`gambit.query.result.QueryInput` in results
		object. Only used for reporting, does not any other aspect of results. Items can be
		``QueryInput``, ``SequenceFile`` or ``str``.
	chunksize
		Number of reference signatures to process at a time. ``None`` means no chunking is performed.
	progress
		Report progress for distance matrix calculation. See
		:func:`gambit.util.progress.get_progress` for description of allowed values.
	"""
	queries = list(queries)
	pconf = progress_config(progress)

	if len(queries) == 0:
		raise ValueError('Must supply at least one query.')

	if inputs is not None:
		inputs = list(map(QueryInput.convert, inputs))
		if len(inputs) != len(queries):
			raise ValueError('Number of inputs does not match number of queries.')
	else:
		inputs = [QueryInput(str(i + 1)) for i in range(len(queries))]

	# Calculate distances
	# (This will only be about 200kB per row/query [50k float32's] so having the whole thing in
	# memory at once isn't a big deal).
	dmat = jaccard_sparse_matrix(
		queries,
		db.signatures,
		ref_indices=db.sig_indices,
		distance=True,
		chunksize=chunksize,
		progress=pconf.update(desc='Calculating distances'),
	)

	with iter_progress(inputs, pconf, desc='Classifying') as inputs_iter:
		items = [classify_item(db, input, dmat[i, :]) for i, input in enumerate(inputs_iter)]

	return QueryResults(items=items, genomeset=db.genomeset, signaturesmeta=db.signatures.meta)


def runquery_parse(db: GAMBITDatabase,
                   files: Sequence[SequenceFile],
                   file_labels: Optional[Sequence[str]] = None,
                   **kw,
                   ) -> QueryResults:
	"""Query a database with signatures derived by parsing a set of genome sequence files.

	Parameters
	----------
	db
		Database to query.
	files
		Sequence files containing query files.
	file_labels
		Custom to use for each file in returned results object. If None will use file names.
	\\**kw
		Additional keyword arguments passed to :func:`.runquery`\\ .
	"""
	from gambit.io.seq import find_kmers_in_files

	pconf = progress_config(kw.pop('progress', None))

	if file_labels is None:
		inputs = files
	else:
		inputs = [QueryInput(label, file) for label, file in zip_strict(file_labels, files)]

	query_sigs = find_kmers_in_files(
		db.signatures.kmerspec,
		files,
		progress=pconf.update(desc='Parsing input'),
	)

	return runquery(db, query_sigs, inputs, progress=pconf, **kw)


def classify_item(db: GAMBITDatabase, input: QueryInput, dists: np.ndarray) -> QueryResultItem:
	"""
	Perform taxonomic classification of single input from distances to reference genomes and return
	a results item object.
	"""
	matches = find_matches(zip_strict(db.genomes, dists))
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
