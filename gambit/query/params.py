from attr import attrs, attrib


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
