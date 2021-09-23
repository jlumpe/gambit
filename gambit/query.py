"""Run queries against a GAMBIT database to predict taxonomy of genome sequences."""
from warnings import warn
from datetime import datetime
from typing import Sequence, Optional, Union, List, Dict, Any

from attr import attrs, attrib

from gambit import __version__ as GAMBIT_VERSION
from gambit.classify import classify, reportable_taxon, ClassifierResult
from gambit.db import GAMBITDatabase, Taxon, ReferenceGenomeSet
from gambit.io.seq import SequenceFile
from gambit.kmers import KmerSignature
from gambit.metric import jaccard_sparse_matrix
from gambit.signatures import SignaturesMeta
from gambit.util.misc import zip_strict
from gambit.util.progress import progress_config, iter_progress


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
	"""
	classify_strict: bool = attrib(default=False)
	chunksize: int = attrib(default=1000)


@attrs()
class QueryInput:
	"""Information on a query genome.

	Attributes
	----------
	label
		Some unique label for the input, probably the file name.
	file
		Source file (optional).
	"""
	label: str = attrib()
	file: Optional[SequenceFile] = attrib(default=None, repr=False)

	@classmethod
	def convert(cls, x: Union['QueryInput', SequenceFile, str]) -> 'QueryInput':
		"""Convenience function to convert flexible argument types into QueryInput.

		Accepts single string label, ``SequenceFile`` (uses file name for label), or existing
		``QueryInput`` instance (returned unchanged).
		"""
		if isinstance(x, QueryInput):
			return x
		if isinstance(x, str):
			return QueryInput(x)
		if isinstance(x, SequenceFile):
			return QueryInput(x.path.name, x)
		raise TypeError(f'Cannot convert {type(x)} instance to QueryInput')


@attrs()
class QueryResultItem:
	"""Result for a single query sequence.

	Attributes
	----------
	input
		Information on input genome.
	classifier_result
		Result of running classifier.
	report_taxon
		Final taxonomy prediction to be reported to the user.
	"""
	input: QueryInput = attrib()
	classifier_result: ClassifierResult = attrib()
	report_taxon: Optional[Taxon] = attrib(default=None)


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
	items: List[QueryResultItem] = attrib()
	params: Optional[QueryParams] = attrib(default=None)
	genomeset: Optional[ReferenceGenomeSet] = attrib(default=None)
	signaturesmeta: Optional[SignaturesMeta] = attrib(default=None)
	gambit_version: str = attrib(default=GAMBIT_VERSION)
	timestamp: datetime = attrib(factory=datetime.now)
	extra: Dict[str, Any] = attrib(factory=dict)


def query(db: GAMBITDatabase,
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


def query_parse(db: GAMBITDatabase,
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
		Additional keyword arguments passed to :func:`.query`\\ .
	"""
	from gambit.signatures.calc import calc_file_signatures

	pconf = progress_config(kw.pop('progress', None))

	if file_labels is None:
		inputs = files
	else:
		inputs = [QueryInput(label, file) for label, file in zip_strict(file_labels, files)]

	query_sigs = calc_file_signatures(
		db.signatures.kmerspec,
		files,
		progress=pconf.update(desc='Parsing input'),
	)

	return query(db, query_sigs, params, inputs=inputs, progress=pconf, **kw)
