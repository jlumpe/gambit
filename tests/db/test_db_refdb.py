"""Test gambit.db.refdb."""

from gambit.db.refdb import ReferenceDatabase, locate_db_files, load_db, load_db_from_dir


def check_loaded_db(db):
	# Just type check for now
	assert isinstance(db, ReferenceDatabase)


def test_locate_db_files(testdb):
	db_file, sigs_file = locate_db_files(testdb.paths.root)
	assert db_file == testdb.paths.ref_genomes
	assert sigs_file == testdb.paths.ref_signatures


def test_load_db(testdb):
	db = load_db(testdb.paths.ref_genomes, testdb.paths.ref_signatures)
	check_loaded_db(db)


def test_load_db_from_dir(testdb):
	db = load_db_from_dir(testdb.paths.root)
	check_loaded_db(db)
