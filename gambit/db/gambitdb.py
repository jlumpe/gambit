from typing import Sequence

from gambit.kmers import KmerSpec
from .models import ReferenceGenomeSet, AnnotatedGenome, genomes_by_id_subset
from gambit.signatures import SignatureArray, SignaturesMeta
from gambit.signatures.base import ReferenceSignatures


class GAMBITDatabase:
	"""Object containing reference genomes, their k-mer signatures, and associated data.

	This is all that is needed at runtime to run queries.

	Attributes
	----------
	kmerspec
		``KmerSpec`` used to calculate k-mer signatures.
	genomeset
		Genome set containing genomes.
	signatures_meta
		Metadata for reference genome signatures.
	genomes
		Reference genomes.
	genome_signatures
		K-mer signatures for each genome.

	Parameters
	----------
	genomeset
	signatures
	"""
	kmerspec: KmerSpec
	genomeset: ReferenceGenomeSet
	signatures_meta: SignaturesMeta
	genomes: Sequence[AnnotatedGenome]
	genome_signatures: SignatureArray

	def __init__(self, genomeset: ReferenceGenomeSet, signatures: ReferenceSignatures):
		self.genomeset = genomeset
		self.kmerspec = signatures.kmerspec
		self.signatures_meta = signatures.meta

		id_attr = signatures.meta.id_attr
		if id_attr is None:
			raise TypeError('id_attr field of signatures metadata cannot be None')

		self.genomes, sig_indices = genomes_by_id_subset(genomeset, id_attr, signatures.ids)
		self.genome_signatures = signatures[sig_indices]

		n = genomeset.genomes.count()
		if len(self.genomes) != n:
			missing = n - len(self.genomes)
			raise ValueError(f'{missing} of {n} genomes not matched to signature IDs. Is the id_attr attribute of the signatures metadata correct?')
