from typing import Sequence

from sqlalchemy.orm import Session, object_session

from .models import ReferenceGenomeSet, AnnotatedGenome, genomes_by_id_subset
from gambit.sigs.meta import ReferenceSignatures


class GAMBITDatabase:
	"""Object containing reference genomes, their k-mer signatures, and associated data.

	This is all that is needed at runtime to run queries.

	Attributes
	----------
	genomeset
		Genome set containing reference genomes.
	genomes
		List of reference genomes.
	signatures
		K-mer signatures for each genome. A subtype of ``ReferenceSignatures``, so contains metadata
		on signatures as well as the signatures themselves. Type may represent signatures stored on
		disk (e.g. :class:`HDF5Signatures`) instead of in memory. OK to contain additional
		signatures not corresponding to any genome in ``genomes``.
	sig_indices
		Index of signature in ``signatures`` corresponding to each genome in ``genomes``.
		In sorted order to improve performance when iterating over them (improve locality if in
		memory and avoid seeking if in file).
	session
		The SQLAlchemy session ``genomeset`` and the elements of ``genomes`` belong to.
		It is important to keep a reference to this, just having references to the ORM objects
		themselves is not enough to keep the session from being garbage collected.

	Parameters
	----------
	genomeset
	signatures
	"""
	genomeset: ReferenceGenomeSet
	genomes: Sequence[AnnotatedGenome]
	signatures: ReferenceSignatures
	sig_indices: Sequence[int]

	def __init__(self, genomeset: ReferenceGenomeSet, signatures: ReferenceSignatures):
		self.genomeset = genomeset
		self.signatures = signatures
		self.session = object_session(genomeset)

		id_attr = signatures.meta.id_attr
		if id_attr is None:
			raise TypeError('id_attr field of signatures metadata cannot be None')

		self.genomes, self.sig_indices = genomes_by_id_subset(genomeset, id_attr, signatures.ids)

		n = genomeset.genomes.count()
		if len(self.genomes) != n:
			missing = n - len(self.genomes)
			raise ValueError(f'{missing} of {n} genomes not matched to signature IDs. Is the id_attr attribute of the signatures metadata correct?')
