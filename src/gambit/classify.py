"""Classify queries based on distance to reference sequences."""

from typing import Optional, Iterable, Sequence

from attr import attrs, attrib
import numpy as np

from gambit.db import AnnotatedGenome, Taxon
from gambit.util.misc import zip_strict


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


def find_matches(itr: Iterable[tuple[AnnotatedGenome, float]]) -> dict[Taxon, list[int]]:
	"""Find taxonomy matches given distances from a query to a set of reference genomes.

	Parameters
	----------
	itr
		Iterable over ``(genome, distance)`` pairs.

	Returns
	-------
	dict[Taxon, list[int]]
		Mapping from taxa to indices of genomes matched to them.
	"""
	matches = dict()

	for i, (g, d) in enumerate(itr):
		match = matching_taxon(g.taxon, d)
		if match is not None:
			matches.setdefault(match, []).append(i)

	return matches


def consensus_taxon(taxa: Iterable[Taxon]) -> tuple[Optional[Taxon], set[Taxon]]:
	"""Take a set of taxa matching a query and find a single consensus taxon for classification.

	If a query matches a given taxon, it is expected that there may be matches to some of that
	taxon's ancestors as well. In this case all matched taxa lie in a single lineage and the most
	specific will be the consensus.

	It may also be possible for a query to match multiple taxa which are "inconsistent" with each
	other in the sense that one is not a descendant of another. In that case the consensus will be
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
	tuple[Optional[Taxon], set[Taxon]]
		Consensus taxon along with the set of any taxa in the argument which are descended from it.
	"""
	taxa = list(taxa)

	# Edge case - input is empty
	if not taxa:
		return (None, set())

	# Current consensus and ancestors, bottom to top
	trunk = list(taxa[0].ancestors(incself=True))

	for taxon in taxa[1:]:
		# Taxon in current trunk, nothing to do
		if taxon in trunk:
			continue

		# Find where ancestry of taxon meets current trunk
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
				# Meets the trunk further up - intersection is new consensus
				trunk = trunk[i:]

			break

		else:
			# No common ancestor exists
			return (None, set(taxa))

	others = {t for t in taxa if t not in trunk}
	return (trunk[0], others)


@attrs()
class GenomeMatch:
	"""Match between a query and a single reference genome.

	This is just used to report the distance from a query to some significant reference genome, it
	does not imply that this distance was close enough to actually make a taxonomy prediction or
	that the prediction was the primary prediction overall.

	Attributes
	----------
	genome
		Reference genome matched to.
	distance
		Distance between query and reference genome.
	matched_taxon
		Taxon prediction based off of this match alone. Will always be ``genome.taxon`` or one of
		its ancestors.
	"""
	genome: AnnotatedGenome = attrib()
	distance: float = attrib()
	matched_taxon: Optional[Taxon] = attrib()

	@matched_taxon.default
	def _matched_taxon_default(self):
		return matching_taxon(self.genome.taxon, self.distance)

	def next_taxon(self) -> Optional[Taxon]:
		"""Get next most specific taxon in lineage of ``genome`` for which the threshold was not met."""
		lo = None
		hi = self.genome.taxon  # Leaf should always have threshold

		while hi is not None:
			if hi.distance_threshold is not None and self.distance <= hi.distance_threshold:
				return lo

			lo = hi

			# Advance to next in ancestry with distance threshold
			hi = hi.parent
			while hi is not None and hi.distance_threshold is None:
				hi = hi.parent

		return lo


@attrs()
class ClassifierResult:
	"""Result of applying the classifier to a single query genome.

	Attributes
	----------
	success
		Whether the classification process ran successfully with no fatal errors. If True it is
		still possible no prediction was made.
	predicted_taxon
		Taxon predicted by classifier.
	primary_match
		Match to closest reference genome which produced a predicted taxon equal to or a descendant
		of ``predicted_taxon``. None if no prediction was made.
	closest_match
		Match to closest reference genome overall. This should almost always be identical to
		``primary_match``.
	next_taxon
		Next most specific taxon for which the threshold was not met. Currently this just taken
		from the ancestry of ``closest_match.genome.taxon``.
	warnings
		List of non-fatal warning messages to report.
	error
		Message describing a fatal error which occurred, if any.
	"""
	success: bool = attrib()
	predicted_taxon: Optional[Taxon] = attrib()
	primary_match: Optional[GenomeMatch] = attrib()
	closest_match: GenomeMatch = attrib()
	next_taxon: Optional[Taxon] = attrib()
	warnings: list[str] = attrib(factory=list, repr=False)
	error: Optional[str] = attrib(default=None, repr=False)

	@next_taxon.default
	def _next_taxon_default(self):
		return self.closest_match.next_taxon()


def classify(ref_genomes: Sequence[AnnotatedGenome],
             dists: np.ndarray,
             *,
             strict: bool = False,
             ) -> ClassifierResult:
	"""Predict the taxonomy of a query genome based on its distances to a set of reference genomes.

	Parameters
	----------
	ref_genomes
		List of reference genomes from database.
	dists
		Array of distances to each reference genome.
	strict
		If true find all significant matches to reference genomes and attempt to reconcile them if
		they result in different taxa. If False just consider the top (closest) match.
		Defaults to False.
	"""
	# Find closest match
	closest = np.argmin(dists)
	closest_match = GenomeMatch(
		genome=ref_genomes[closest],
		distance=dists[closest],
		matched_taxon=matching_taxon(ref_genomes[closest].taxon, dists[closest]),
	)

	if not strict:
		# Use closest match only
		return ClassifierResult(
			success=True,
			predicted_taxon=closest_match.matched_taxon,
			primary_match=closest_match if closest_match.matched_taxon is not None else None,
			closest_match=closest_match,
		)

	# Find all matches and attempt to get consensus
	matches = find_matches(zip_strict(ref_genomes, dists))
	consensus, others = consensus_taxon(matches.keys())

	# No matches found
	if not matches:
		return ClassifierResult(
			success=True,
			predicted_taxon=None,
			primary_match=None,
			closest_match=closest_match,
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
			genome=ref_genomes[best_i],
			distance=best_d,
			matched_taxon=best_taxon,
		)

	result = ClassifierResult(
		success=True,
		predicted_taxon=consensus,
		primary_match=primary_match,
		closest_match=closest_match,
	)

	# Warn of inconsistent matches
	if others:
		msg = f'Query matched {len(others)} inconsistent taxa: '
		msg += ', '.join(sorted(other.short_repr() for other in others))
		msg += '. Reporting lowest common ancestor of this set.'
		result.warnings.append(msg)

	# No consensus found - matches do not have common ancestor
	if consensus is None:
		result.success = False
		result.error = 'Matched taxa have no common ancestor.'

	# Primary match is not closest
	if primary_match is not None and primary_match.genome != closest_match.genome:
		result.warnings.append('Primary genome match is not closest match.')

	return result
