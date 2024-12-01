"""Run queries against a GAMBIT database to predict taxonomy of genome sequences."""

from warnings import warn
from datetime import datetime
from typing import Sequence, Optional, Any
from pathlib import Path

from attr import attrs, attrib
from attr.converters import optional as optional_converter
import numpy as np

from gambit import __version__ as GAMBIT_VERSION
from gambit.classify import classify, ClassifierResult, GenomeMatch
from gambit.db import ReferenceDatabase, Taxon, ReferenceGenomeSet, reportable_taxon
from gambit.sigs.base import KmerSignature, SignaturesMeta, ReferenceSignatures
from gambit.metric import jaccarddist_matrix
from gambit.util.io import FilePath
from gambit.util.progress import progress_config, iter_progress
from gambit.util.misc import zip_strict


@attrs()
class QueryParams:
	"""Parameters for running a query.

	Attributes
	----------
	classify_strict
		``strict`` parameter to :func:`gambit.classify.classify`. Defaults to False.
	chunksize
		Number of reference signatures to process at a time. ``None`` means no chunking is performed.
		Defaults to 1000.
	report_closest
		Number of closest genomes to report in results. Does not affect classification.
	"""
	classify_strict: bool = attrib(default=False)
	chunksize: int = attrib(default=1000)
	report_closest: int = attrib(default=10)


@attrs()
class QueryResultItem:
	"""Result for a single query sequence.

	Attributes
	----------
	label
		Unique label describing query.
	classifier_result
		Result of running classifier.
	report_taxon
		Final taxonomy prediction to be reported to the user.
	closest_genomes
		List of closest reference genomes to query. Length determined by
		:attr:`.QueryParams.report_closest`.
	file
		Path to file containing query genome (optional).
	"""
	label: str = attrib()
	classifier_result: ClassifierResult = attrib()
	report_taxon: Optional[Taxon] = attrib(default=None)
	closest_genomes: list[GenomeMatch] = attrib(factory=list)
	file: Optional[Path] = attrib(default=None, converter=optional_converter(Path))


@attrs(repr=False)
class QueryResults:
	"""Results for a set of queries, as well as information on database and parameters used.

	Attributes
	----------
	items
		Results for each query sequence.
	params
		Parameters used to run query.
	genomeset
		Genome set used.
	signaturesmeta
		Metadata for signatures set used.
	gambit_version
		Version of GAMBIT command/library used to generate the results.
	timestamp
		Time query was completed.
	extra
		JSON-able dict containing additional arbitrary metadata.
	"""
	items: list[QueryResultItem] = attrib()
	params: Optional[QueryParams] = attrib(default=None)
	genomeset: Optional[ReferenceGenomeSet] = attrib(default=None)
	signaturesmeta: Optional[SignaturesMeta] = attrib(default=None)
	gambit_version: str = attrib(default=GAMBIT_VERSION)
	timestamp: datetime = attrib(factory=datetime.now)
	extra: dict[str, Any] = attrib(factory=dict)


def query(db: ReferenceDatabase,
          queries: Sequence[KmerSignature],
          params: Optional[QueryParams] = None,
          *,
          labels: Optional[Sequence[str]] = None,
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
		``QueryParams`` instance defining parameter values. If None take values from additional
		keyword arguments or use defaults.
	labels
		Optional list of string labels for each query. Only used for reporting (sets ``label``
		attribute of :class:`QueryResultItem` in results object), does not any other aspect of
		results.
	progress
		Report progress for distance matrix calculation and classification. See
		:func:`gambit.util.progress.get_progress` for description of allowed values.
	\\**kw
		Passed to ``QueryParams``.
	"""
	if params is None:
		params = QueryParams(**kw)
	elif kw:
		warn('Additional keyword arguments ignored if "params" argument is not None.')

	pconf = progress_config(progress)

	if len(queries) == 0:
		raise ValueError('Must supply at least one query.')

	# Labels
	if labels is not None:
		if len(labels) != len(queries):
			raise ValueError('Number of labels does not match number of queries.')

	elif isinstance(queries, ReferenceSignatures):
		# Get default labels from queries of ReferenceSignatures object
		labels = list(map(str, queries.ids))

	else:
		labels = [str(i + 1) for i in range(len(queries))]

	# Calculate distances
	# (This will only be about 200kB per row/query [50k float32's] so having the whole thing in
	# memory at once isn't a big deal).
	dmat = jaccarddist_matrix(
		queries,
		db.signatures,
		ref_indices=db.sig_indices,
		chunksize=params.chunksize,
		progress=pconf.update(desc='Calculating distances'),
	)

	# Classify inputs and create result items
	with iter_progress(labels, pconf, desc='Classifying') as labels_iter:
		items = [
			get_result_item(db, params, dmat[i, :], label)
			for i, label in enumerate(labels_iter)
		]

	return QueryResults(
		items=items,
		params=params,
		genomeset=db.genomeset,
		signaturesmeta=db.signatures.meta,
	)


def get_result_item(db: ReferenceDatabase, params: QueryParams, dists: np.ndarray, label: str) -> QueryResultItem:
	"""Perform classification and create result item object for single query input.

	Parameters
	----------
	db
	params
	dists
		1D array of distances from query to all reference genomes.
	label
	"""
	clsresult = classify(db.genomes, dists, strict=params.classify_strict)
	closest = [GenomeMatch(db.genomes[i], dists[i]) for i in np.argsort(dists)[:params.report_closest]]

	return QueryResultItem(
		label=label,
		classifier_result=clsresult,
		report_taxon=reportable_taxon(clsresult.predicted_taxon),
		closest_genomes=closest,
	)


def query_parse(db: ReferenceDatabase,
                files: Sequence[FilePath],
                params: Optional[QueryParams] = None,
                *,
                labels: Optional[Sequence[str]] = None,
                parse_kw: Optional[dict[str, Any]] = None,
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
		``QueryParams`` instance defining parameter values. If None take values from additional
		keyword arguments or use defaults.
	labels
		Custom labels to use for each file in returned results object. If None use file names.
	parse_kw
		Keyword parameters to pass to :func:`gambit.sigs.calc.calc_file_signatures`.
	\\**kw
		Additional keyword arguments passed to :func:`.query`.
	"""
	from gambit.sigs.calc import calc_file_signatures

	pconf = progress_config(kw.pop('progress', None))
	if parse_kw is None:
		parse_kw = dict()
	parse_kw.setdefault('progress', pconf.update(desc='Parsing input'))

	if labels is None:
		labels = [str(file) for file in files]
	else:
		if len(labels) != len(files):
			raise ValueError('Number of labels does not match number of files')

	query_sigs = calc_file_signatures(db.signatures.kmerspec, files, **parse_kw)

	results = query(db, query_sigs, params, labels=labels, progress=pconf, **kw)

	# Assign file attribute of QueryResultItem's
	for item, file in zip_strict(results.items, files):
		item.file = file

	return results
