"""Data classes to represent query results."""

from typing import Optional, List, Union, Dict, Any
from datetime import datetime

from attr import attrs, attrib

from gambit.io.seq import SequenceFile
from gambit.db.models import ReferenceGenomeSet, Taxon
from gambit.signatures import SignaturesMeta
from gambit.classify import ClassifierResult
from .params import QueryParams
from gambit import __version__ as GAMBIT_VERSION


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
