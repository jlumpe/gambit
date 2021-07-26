"""Classify queries based on distance to reference sequences."""

from typing import Optional, Tuple, Iterable, Dict, List, Set

from gambit.db.models import AnnotatedGenome, Taxon


def matching_taxon(taxon: Taxon, d: float) -> Optional[Taxon]:
	"""Find first taxon in linage for which distance ``d`` is within its classification threshold.

	Parameters
	----------
	taxon
		Taxon to start searching from.
	d
		Distance value.

	Returns
	-------
	Optional[Taxon]
		Most specific taxon in ancestry with ``threshold_distance >= d``.
	"""
	for t in taxon.ancestors(incself=True):
		if t.distance_threshold is not None and d <= t.distance_threshold:
			return t
	return None


def find_matches(itr: Iterable[Tuple[AnnotatedGenome, float]]) -> Dict[Taxon, List[int]]:
	"""Find taxonomy matches given distances from a query to a set of reference genomes.

	Parameters
	----------
	itr
		Iterable over ``(genome, distance)`` pairs.

	Returns
	-------
	Dict[Taxon, List[Int]]
		Mapping from taxa to indices of genomes matched to them.
	"""
	matches = dict()

	for i, (g, d) in enumerate(itr):
		match = matching_taxon(g.taxon, d)
		if match is not None:
			matches.setdefault(match, []).append(i)

	return matches


def consensus_taxon(taxa: Iterable[Taxon]) -> Tuple[Optional[Taxon], Set[Taxon]]:
	"""Take a set of taxa matching a query and find a single consensus taxon for classification.

	If a query matches a given taxon, it is expected that there may be hits on some of its ancestors
	as well. In this case all taxa lie in a single lineage and the most specific taxon will be the
	consensus.

	It may also be possible for a query to match multiple taxa which are "inconsistent" with each
	other in the sense that one is not a descendant of the other. In that case the consensus will be
	the lowest taxon which is either a descendant or ancestor of all taxa in the argument. It's also
	possible in pathological cases (depending on reference database design) that the taxa may be
	within entirely different trees, in which case the consensus will be ``None``. The second
	element of the returned tuple is the set of taxa in the argument which are strict descendants of
	the consensus. This set will contain at least two taxa in the case of such an inconsistency and
	be empty otherwise.

	Parameters
	----------
	taxa : Iterable[Taxon]

	Returns
	-------
	Tuple[Optional[Taxon], List[Taxon]]
		Consensus taxon along with subset of ``taxa`` which are not an ancestor of it.
	"""
	taxa = list(taxa)

	# Edge case - input is empty
	if not taxa:
		return (None, set())

	# Current consensus and ancestors, bottom to top
	trunk = list(taxa[0].ancestors(incself=True))
	# Taxa which are strict descendants of current consensus
	crown = set()

	for taxon in taxa[1:]:
		# Taxon in current trunk, nothing to do
		if taxon in trunk:
			continue

		# Find lowest ancestor of taxon in current trunk
		for a in taxon.ancestors(incself=False):
			try:
				i = trunk.index(a)
			except ValueError:
				# Current ancestor not in trunk, continue to parent
				continue

			if i == 0:
				# Directly descended from current consensus, this taxon becomes new consensus
				trunk = list(taxon.ancestors(incself=True))
			else:
				# Descended from one of current consensus's ancestors

				# This taxon and current consensus added to crown
				crown.add(trunk[0])
				crown.add(taxon)

				# New consensus is this ancestor
				trunk = trunk[i:]

			break

		else:
			# No common ancestor exists
			return (None, set(taxa))

	return (trunk[0] if trunk else None, crown)


def reportable_taxon(taxon: Taxon) -> Optional[Taxon]:
	"""Find the first reportable taxon in a linage.

	Parameters
	----------
	taxon
		Taxon to start looking from.

	Returns
	-------
	Optional[gambit.db.models.Taxon]
		Most specific taxon in ancestry with ``report=True``, or ``None`` if none found.
	"""
	for t in taxon.ancestors(incself=True):
		if t.report:
			return t

	return None
