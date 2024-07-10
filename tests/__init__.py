"""
Adding the __init__.py to the tests/ directory (and its subdirectories) makes them all part of the
same package structure.

- Allows test modules/files to import from each other (including from modules in different
  directories, such as files in tests/cli/ importing from tests/testdb.py).
- Does not require test modules to have unique names.

This necessitates using the "prepend" (or possibly "append"?) import mode (which is the default).
This setup comes with its own set of caveats. See
https://docs.pytest.org/en/7.1.x/explanation/pythonpath.html for a discussion of how test modules
are imported.
"""
