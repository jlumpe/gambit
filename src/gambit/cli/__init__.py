"""GAMBIT command line interface."""

from .root import cli as cli

# Imported for side effect of registering commands on .root.cli
from . import query  # noqa: F401
from . import signatures  # noqa: F401
from . import dist  # noqa: F401
from . import tree  # noqa: F401
from . import debug  # noqa: F401
