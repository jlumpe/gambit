"""Test gambit.db.refdb."""

from gambit.db.refdb import ReferenceDatabase, locate_db_files, load_db, load_db_from_dir


def check_loaded_db(db):
	# Just type check for now
	assert isinstance(db, ReferenceDatabase)


def test_locate_db_files(testdb_files):
	db_file, sigs_file = locate_db_files(testdb_files['root'])
	assert db_file == testdb_files['ref_genomes']
	assert sigs_file == testdb_files['ref_signatures']


def test_load_db(testdb_files):
	db = load_db(testdb_files['ref_genomes'], testdb_files['ref_signatures'])
	check_loaded_db(db)


def test_load_db_from_dir(testdb_files):
	db = load_db_from_dir(testdb_files['root'])
	check_loaded_db(db)
