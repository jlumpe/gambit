"""Development tools."""

from pathlib import Path
import subprocess
import shutil
from typing import Dict, Any

import gambit
from gambit.io.util import FilePath
from gambit.util.misc import zip_strict


_INSTALL_INFO = None


def get_commit_info(repo_path: FilePath, commit: str = 'HEAD') -> Dict[str, str]:
	"""Get metadata on a git commit.

	This calls the ``git`` command, so it must be installed and available.

	Parameters
	----------
	repo_path
		Path to git repo.
	commit
		Commit to get information on.
	"""
	fields = [
		('hash', '%H'),
		('author', '%an <%ae>'),
		('author_date', '%aI'),
		('commit', '%cn <%ce>'),
		('commit_date', '%cI'),
		('subject', '%s'),
	]

	fmt_str = '%n'.join(fmt for name, fmt in fields)
	cmd = ['git', 'show', '-s', '--format=' + fmt_str, commit]

	result = subprocess.run(cmd, cwd=repo_path, capture_output=True, check=True, text=True)

	lines = result.stdout.splitlines()
	assert len(lines) == len(fields)
	return {name: line for (name, fmt), line in zip_strict(fields, lines)}


def _install_info():
	info = dict(pkg_dir=None, repo_dir=None, commit=None)

	if not hasattr(gambit, '__path__'):
		info['status'] = 'gambit module has no __path__ attribute.'
		return info

	if len(gambit.__path__) != 1:
		info['status'] = f'Expected gambit.__path__ to contain single item, got {gambit.__path__!r}'
		return info

	pkg_dir = info['pkg_dir'] = Path(gambit.__path__[0])
	repo_dir = pkg_dir.parent

	if not (repo_dir / '.git').is_dir():
		info['status'] = 'Parent of package directory not a git repo (has no .git subdirectory).'
		return info

	info['repo_dir'] = repo_dir

	if shutil.which('git') is None:
		info['status'] = 'git command not found'
		return info

	try:
		commit = get_commit_info(repo_dir)
	except subprocess.SubprocessError as e:
		info['status'] = f'Command {e.cmd!r} returned exit code {e.returncode} with stderr output {e.stderr!r}'
	except Exception as e:
		info['status'] = f'Error getting commit info: {e!r}'
	else:
		info['status'] = 'Git info retrieved successfully.'
		info['commit'] = commit

	return info


def install_info() -> Dict[str, Any]:
	"""Get information on the GAMBIT installation if it is installed in development mode.


	If gambit is installed via the setuptools development install method (``pip install -e``), this
	checks if the source directory is a valid git repo and tries to get information on the current
	commit. This is used to mark exported results from development versions of the software which do
	not correspond to an official release.
	"""
	global _INSTALL_INFO
	if _INSTALL_INFO is None:
		_INSTALL_INFO = _install_info()
	return _INSTALL_INFO
