"""Query a GAMBIT database."""

from typing import Sequence, Optional, Union
from warnings import warn

from gambit.db.gambitdb import GAMBITDatabase
from gambit.kmers import KmerSignature
from gambit.io.seq import SequenceFile
from gambit.util.misc import zip_strict
from gambit.util.progress import progress_config, iter_progress
from gambit.metric import jaccard_sparse_matrix
from gambit.classify import classify, reportable_taxon
from .params import QueryParams
from .results import QueryInput, QueryResultItem, QueryResults


def runquery(db: GAMBITDatabase,
             queries: Sequence[KmerSignature],
             params: Optional[QueryParams] = None,
             *,
             inputs: Optional[Sequence[Union[QueryInput, SequenceFile, str]]] = None,
             progress = None,
             **kw,
             ) -> QueryResults:
	"""Predict the taxonomy of one or more query genomes using a GAMBIT reference database.

	Parameters
	----------
	db
		Database to query.
	queries
		Sequence of k-mer signatures for query genomes.
	params
		``QueryParams`` instance defining parameter values. If None will take values from additional
		keyword arguments or use defaults.
	inputs
		Description for each input, converted to :class:`gambit.query.result.QueryInput` in results
		object. Only used for reporting, does not any other aspect of results. Items can be
		``QueryInput``, ``SequenceFile`` or ``str``.
	progress
		Report progress for distance matrix calculation. See
		:func:`gambit.util.progress.get_progress` for description of allowed values.
	\\**kw
		Passed to ``QueryParams``\\ .
	"""
	if params is None:
		params = QueryParams(**kw)
	elif kw:
		warn('Additional keyword arguments ignored if "params" argument is not None.')

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
		chunksize=params.chunksize,
		progress=pconf.update(desc='Calculating distances'),
	)

	items = []

	# Classify inputs and create result items
	with iter_progress(inputs, pconf, desc='Classifying') as inputs_iter:
		for i, input in enumerate(inputs_iter):
			clsresult = classify(db.genomes, dmat[i, :], strict=params.classify_strict)
			report_taxon = None if clsresult.predicted_taxon is None else reportable_taxon(clsresult.predicted_taxon)
			items.append(QueryResultItem(
				input=input,
				classifier_result=clsresult,
				report_taxon=report_taxon,
			))

	return QueryResults(
		items=items,
		params=params,
		genomeset=db.genomeset,
		signaturesmeta=db.signatures.meta,
	)


def runquery_parse(db: GAMBITDatabase,
                   files: Sequence[SequenceFile],
                   params: Optional[QueryParams] = None,
                   *,
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
	params
		``QueryParams`` instance defining parameter values. If None will take values from additional
		keyword arguments or use defaults.
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

	return runquery(db, query_sigs, params, inputs=inputs, progress=pconf, **kw)
