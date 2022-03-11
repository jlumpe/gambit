"""Test gambit.db.refdb."""

from gambit.db.refdb import ReferenceDatabase


def check_loaded_db(db):
	# Just type check for now
	assert isinstance(db, ReferenceDatabase)


class TestReferenceDatabase:
	"""Test the ReferenceDatabase class."""

	def test_locate_files(self, testdb):
		db_file, sigs_file = ReferenceDatabase.locate_files(testdb.paths.root)
		assert db_file == testdb.paths.ref_genomes
		assert sigs_file == testdb.paths.ref_signatures

	def test_load(self, testdb):
		db = ReferenceDatabase.load(testdb.paths.ref_genomes, testdb.paths.ref_signatures)
		check_loaded_db(db)

	def test_load_db_from_dir(self, testdb):
		db = ReferenceDatabase.load_from_dir(testdb.paths.root)
		check_loaded_db(db)
