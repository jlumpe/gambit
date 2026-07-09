"""
Genomic Approximation Method for Bacterial Identification and Tracking.

Author: Jared Lumpe
Author email: jared@jaredlumpe.com
"""

__author__ = 'Jared Lumpe'


from typing import Any


# Previously the package version was defined here in the __version__ variable. It has since been
# moved to pyproject.toml. Provide this for backwards compatibility.
def __getattr__(name: str) -> Any:
	if name == '__version__':
		from importlib.metadata import version, PackageNotFoundError
		try:
			return version('gambit')
		except PackageNotFoundError:
			return None
	raise AttributeError(f"module '{__name__}' has no attribute '{name}'")
