"""Test gambit.db.sqla."""

from sqlalchemy.orm import Session

from gambit.db import ReadOnlySession, file_sessionmaker
from ..testdb import TestDB


def test_file_sessionmaker(testdb: TestDB):
	db_file = testdb.paths.ref_genomes

	maker = file_sessionmaker(db_file, readonly=True)
	assert isinstance(maker(), ReadOnlySession)

	maker = file_sessionmaker(db_file, readonly=False)
	assert isinstance(maker(), Session)

	for cls in [Session, ReadOnlySession]:
		maker = file_sessionmaker(db_file, class_=cls)
		assert isinstance(maker(), cls)
