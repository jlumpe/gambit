"""Query a GAMBIT database."""

from typing import Sequence, Optional, Union

import numpy as np

from gambit.db.gambitdb import GAMBITDatabase
from gambit.kmers import KmerSignature
from gambit.io.seq import SequenceFile
from gambit.util.misc import zip_strict
from gambit.util.progress import progress_config, iter_progress
from gambit.metric import jaccard_sparse_matrix
from gambit.classify import classify, reportable_taxon
from .results import QueryInput, QueryResultItem, QueryResults


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

	items = []

	# Classify inputs and create result items
	with iter_progress(inputs, pconf, desc='Classifying') as inputs_iter:
		for i, input in enumerate(inputs_iter):
			clsresult = classify(db.genomes, dmat[i, :])
			report_taxon = None if clsresult.predicted_taxon is None else reportable_taxon(clsresult.predicted_taxon)
			items.append(QueryResultItem(
				input=input,
				classifier_result=clsresult,
				report_taxon=report_taxon,
			))

	return QueryResults(
		items=items,
		genomeset=db.genomeset,
		signaturesmeta=db.signatures.meta,
	)


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
