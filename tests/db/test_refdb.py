"""Test gambit.db.refdb."""

import random

import pytest
from sqlalchemy.orm import sessionmaker

from gambit.db import refdb
from gambit.db import Genome, ReferenceGenomeSet, AnnotatedGenome, Taxon, ReferenceDatabase, DatabaseLoadError


GENOME_ID_ATTRS = {attr: getattr(Genome, attr) for attr in Genome.ID_ATTRS}


def check_loaded_db(db):
	# Just type check for now
	assert isinstance(db, ReferenceDatabase)


class TestGenomeIDMapping:
	"""Test mapping genomes to ID values."""

	@pytest.fixture()
	def session(self, make_empty_db):
		"""In-memory database containing genomes which have values for all ID attributes."""
		engine = make_empty_db()
		Session = sessionmaker(engine)
		session = Session()

		gset = ReferenceGenomeSet(
			key='test_gset',
			version='1.0',
			name='Test genome set',
		)
		session.add(gset)

		roottaxon = Taxon(
			key='root',
			name='root',
			genome_set=gset,
		)
		session.add(roottaxon)

		for i in range(20):
			g = Genome(
				key=f'genome_{i}',
				description=f'Test genome {i}',
				ncbi_db='assembly',
				ncbi_id=i,
				genbank_acc=f'GCA_{i:09d}.1',
				refseq_acc=f'GCF_{i:09d}.1',
			)
			ag = AnnotatedGenome(
				genome_set=gset,
				genome=g,
				taxon=roottaxon,
			)
			session.add(ag)

		session.commit()
		return session

	def test__genome_id_attr(self):
		"""Test _check_genome_id_attr() function."""

		for key, attr in GENOME_ID_ATTRS.items():
			assert refdb._check_genome_id_attr(key) is attr
			assert refdb._check_genome_id_attr(attr) is attr

		for arg in ['description', Genome.description, AnnotatedGenome.key]:
			with pytest.raises(ValueError):
				refdb._check_genome_id_attr(arg)

	def test__get_genome_id(self, session):
		"""Test _get_genome_id() function."""

		for genome in session.query(Genome):
			for key, attr in GENOME_ID_ATTRS.items():
				assert refdb._get_genome_id(genome, attr) == getattr(genome, key)

	def test_genomes_by_id(self, session):
		"""Test genomes_by_id() and genomes_by_id_subset() functions."""
		random.seed(0)

		gset = session.query(ReferenceGenomeSet).one()

		# Desired result with only valid IDs
		genomes = list(gset.genomes)
		random.shuffle(genomes)

		# Desired result with some invalid IDs mixed in
		genomes_missing = list(genomes)
		genomes_missing.insert(3, None)
		genomes_missing.insert(8, None)

		for key, attr in GENOME_ID_ATTRS.items():
			# Check with all valid IDs
			ids = [refdb._get_genome_id(g, attr) for g in genomes]
			assert refdb.genomes_by_id(gset, attr, ids) == genomes
			assert refdb.genomes_by_id(gset, key, ids) == genomes

			# Add IDs which cannot be mapped to genomes
			ids_missing = [
				refdb._get_genome_id(g, attr) if g is not None else '!INVALID!'
				for g in genomes_missing
			]

			# Raises KeyError with strict=True (default)
			with pytest.raises(KeyError):
				refdb.genomes_by_id(gset, attr, ids_missing)

			# strict=False returns list with None values
			assert refdb.genomes_by_id(gset, attr, ids_missing, strict=False) == genomes_missing

			# Test genomes_by_id_subset
			genomes_sub, idxs_sub = refdb.genomes_by_id_subset(gset, attr, ids_missing)
			assert genomes_sub == genomes
			assert [genomes_missing[i] for i in idxs_sub] == genomes

			# Incomplete set of IDs which does not encompass all genomes
			ids_incomplete = ids[:-1]

		# TODO


class TestReferenceDatabase:
	"""Test the ReferenceDatabase class."""

	def test_locate_files(self, tmp_path):
		genomes = tmp_path / 'test.gdb'
		genomes2 = tmp_path / 'test2.gdb'
		signatures = tmp_path / 'test.gs'

		# None
		with pytest.raises(DatabaseLoadError):
			ReferenceDatabase.locate_files(tmp_path)

		# Genomes but no signatures
		genomes.touch()
		with pytest.raises(DatabaseLoadError):
			ReferenceDatabase.locate_files(tmp_path)

		# Both
		signatures.touch()
		assert ReferenceDatabase.locate_files(tmp_path) == (genomes, signatures)

		# Extra genomes file
		genomes2.touch()
		with pytest.raises(DatabaseLoadError):
			ReferenceDatabase.locate_files(tmp_path)

		# Alternate extensions
		for file in [genomes, genomes2, signatures]:
			file.unlink()
		genomes = tmp_path / 'test.db'
		signatures = tmp_path / 'test.gs'
		genomes.touch()
		signatures.touch()
		assert ReferenceDatabase.locate_files(tmp_path) == (genomes, signatures)

	def test_load(self, testdb):
		db = ReferenceDatabase.load(testdb.paths.ref_genomes, testdb.paths.ref_signatures)
		check_loaded_db(db)

	def test_load_db_from_dir(self, testdb):
		db = ReferenceDatabase.load_from_dir(testdb.paths.root)
		check_loaded_db(db)
