import os
from typing import Optional, Sequence, TextIO, Union, Iterable, Tuple, List
from pathlib import Path
from collections import Counter

import click
from sqlalchemy import create_engine
from sqlalchemy.engine import Engine
from sqlalchemy.orm import sessionmaker

from gambit.kmers import KmerSpec, DEFAULT_KMERSPEC
from gambit.db import ReferenceDatabase, ReadOnlySession, only_genomeset, DatabaseLoadError
from gambit.sigs.base import ReferenceSignatures, load_signatures
from gambit.util.io import FilePath, read_lines
from gambit.util.misc import join_list_human
from gambit.seq import validate_dna_seq_bytes, SequenceFile


class CLIContext:
	"""Click context object for GAMBIT CLI.

	Loads reference database data lazily the first time it is requested.

	Currently a single option (or environment variable) is used to specify the location of the
	database files, in the future options may be added to specify the reference genomes SQLite
	file and genome signatures file separately. Class methods treat them as being independent.

	Attributes
	----------
	root_context
		Click context object from root command group.
	db_path
		Path to directory containing database files, specified in root command group.
	has_genomes
		Whether reference genome metadata is available.
	has_signatures
		Whether reference signatures are available.
	has_database
		Whether reference genome metadata and reference signatures are both available.
	engine
		SQLAlchemy engine connecting to genomes database.
	Session
		SQLAlchemy session maker for genomes database.
	signatures
		Reference genome signatures.
	"""
	root_context: click.Context
	db_path: Optional[Path]
	has_genomes: bool
	has_signatures: bool
	has_database: bool
	engine: Optional[Engine]
	Session: Optional[sessionmaker]
	signatures: Optional[ReferenceSignatures]

	def __init__(self, root_context: click.Context):
		"""
		Parameters
		----------
		root_context
			Click context object from root command group.
		"""
		self.root_context = root_context

		db_path = root_context.params['db_path']
		self.db_path = None if db_path is None else Path(db_path)

		self._db_found = False
		self._has_genomes = None
		self._has_signatures = None
		self._signatures_path = None

		self._engine = None
		self._Session = None
		self._signatures = None

	def _find_db(self):
		"""Find database files."""
		if self._db_found:
			return

		if self.db_path is None:
			self._has_genomes = self._has_signatures = False

		else:
			try:
				self._genomes_path, self._signatures_path = ReferenceDatabase.locate_files(self.db_path)
			except DatabaseLoadError as e:
				raise click.ClickException(str(e))
			self._has_genomes = self._has_signatures = True

		self._db_found = True

	@property
	def has_genomes(self):
		if not self._db_found:
			self._find_db()
		return self._has_genomes

	@property
	def has_signatures(self):
		if not self._db_found:
			self._find_db()
		return self._has_signatures

	@property
	def has_database(self):
		return self.has_genomes and self.has_signatures

	def require_database(self):
		"""Raise an exception if genome metadata and signatures are not available."""
		if not self.has_database:
			raise click.ClickException('Must supply path to database directory.')

	def require_genomes(self):
		"""Raise an exception if genome metadata is not available."""
		self.require_database()

	def require_signatures(self):
		"""Raise an exception if signatures are not available."""
		self.require_database()

	def _init_genomes(self):
		if self._engine is not None or not self.has_genomes:
			return

		self._engine = create_engine(f'sqlite:///{self._genomes_path}')
		self._Session = sessionmaker(self.engine, class_=ReadOnlySession)

	@property
	def engine(self):
		self._init_genomes()
		return self._engine

	@property
	def Session(self):
		self._init_genomes()
		return self._Session

	@property
	def signatures(self):
		if self._signatures is None and self.has_signatures:
			self._signatures = load_signatures(self._signatures_path)

		return self._signatures

	def get_database(self) -> ReferenceDatabase:
		"""Get reference database object."""
		self.require_database()
		session = self.Session()
		gset = only_genomeset(session)
		return ReferenceDatabase(gset, self.signatures)


################################################################################
# Shared CLI parameters
################################################################################

def filepath(**kw) -> click.Path:
	"""Click Path argument type accepting files only."""
	kw.setdefault('path_type', Path)
	return click.Path(file_okay=True, dir_okay=False, **kw)


def dirpath(**kw) -> click.Path:
	"""Click Path argument type accepting directories only."""
	kw.setdefault('path_type', Path)
	return click.Path(file_okay=False, dir_okay=True, **kw)


def genome_files_arg():
	"""Click positional argument for genome files."""
	return click.argument(
		'files_arg',
		nargs=-1,
		type=filepath(exists=True),
		metavar='GENOMES...',
	)


def cores_param():
	"""Click parameter for number of CPU cores."""
	return click.option('-c', '--cores', type=click.IntRange(min=1), help='Number of CPU cores to use.')


def progress_param():
	"""Click argument to show progress meter."""
	return click.option('--progress/--no-progress', default=True, help="Show/don't show progress meter.")


def listfile_param(*param: str, **kw):
	"""Returns decorator to add param for file listing input paths."""
	return click.option(*param, type=click.File('r'), **kw)


def listfile_dir_param(*param: str, file_metavar=None, **kw):
	"""Returns decorator to add param for parent directory of paths in list file."""
	kw.setdefault('default', '.')
	if file_metavar is not None:
		kw.setdefault('help', f'Parent directory of paths in {file_metavar}.')

	return click.option(*param, type=dirpath(), **kw)


def kspec_params(default: bool = False):
	"""Returns a decorator to add k and prefix options to command.

	Parameters
	----------
	default
		Whether to add default values.
	"""
	popt = click.option(
		'-p', '--prefix',
		help='K-mer prefix.',
		default=DEFAULT_KMERSPEC.prefix_str if default else None,
		metavar='NUCS',
	)
	kopt = click.option(
		'-k',
		type=int,
		help='Number of nucleotides to recognize AFTER prefix.',
		default=DEFAULT_KMERSPEC.k if default else None,
	)
	return lambda f: kopt(popt(f))


def kspec_from_params(k: Optional[int], prefix: Optional[str], default: bool = False) -> Optional[KmerSpec]:
	"""Get KmerSpec from CLI arguments and validate.

	Parameters
	----------
	k
	prefix
	default
		Return default KmerSpec if arguments are None.
	"""

	if prefix is None and k is None:
		return DEFAULT_KMERSPEC if default else None

	if prefix is None or k is None:
		raise click.ClickException('Must specify values for both -k and --prefix arguments.')

	# TODO - minimum k and prefix length are fairly arbitrary here - is there a better method?
	MIN_K = 5
	if k < MIN_K:
		raise click.ClickException(f'k must be at least {MIN_K}')

	MIN_PREFIX_LEN = 2
	if len(prefix) < MIN_PREFIX_LEN:
		raise click.ClickException(f'Prefix length must be at least {MIN_PREFIX_LEN}')

	prefix_bytes = prefix.upper().encode('ascii')
	try:
		validate_dna_seq_bytes(prefix_bytes)
	except ValueError:
		raise click.ClickException(f'Invalid nucleotide codes in prefix: {prefix}')

	return KmerSpec(k, prefix_bytes)


################################################################################
# Sequence file input
################################################################################

FASTA_EXTENSIONS = ('.fasta', '.fna', '.ffn', '.faa', '.frn', '.fa')
GZIP_EXTENSIONS = ('.gz',)


def strip_extensions(filename: str, extensions: Iterable[str]) -> str:
	for ext in extensions:
		if filename.endswith(ext):
			return filename[:-len(ext)]
	return filename


def strip_seq_file_ext(filename: str) -> str:
	"""Strip FASTA and/or gzip extensions from sequence file name."""
	filename = strip_extensions(filename, GZIP_EXTENSIONS)
	filename = strip_extensions(filename, FASTA_EXTENSIONS)
	return filename


def get_file_id(path: FilePath, strip_dir: bool = True, strip_ext: bool = True) -> str:
	"""Get sequence file ID derived from file path.

	Parameters
	----------
	strip_dir
		Strip leading path components.
	strip_ext
		Strip file extension(s).
	"""
	id = os.fspath(path)
	if strip_dir:
		id = os.path.basename(id)
		if strip_ext:
			id = strip_seq_file_ext(id)
	return id


def get_sequence_files(explicit: Optional[Iterable[FilePath]]=None,
                       listfile: Union[None, FilePath, TextIO]=None,
                       listfile_dir: Optional[str]=None,
                       strip_dir: bool = True,
                       strip_ext: bool = True,
                       ) -> Union[Tuple[List[str], List[SequenceFile]], Tuple[None, None]]:
	"""Get list of sequence file paths and IDs from several types of CLI arguments.

	Does not check for conflict between ``explicit`` and ``listfile``.

	Parameters
	----------
	explicit
		List of paths given explicitly, such as with a positional argument.
	listfile
		File listing sequence files, one per line.
	listfile_dir
		Parent directory for files in ``listfile``.
	strip_dir
		Strip leading path components from file paths to derive IDs.
	strip_ext
		Strip file extension from file names to derive IDs.

	Returns
	-------
	Tuple[Optional[List[str]], Optional[List[SequenceFile]]]
		``(ids, files)`` tuple. ``ids`` is a list of string IDs that can be used to label output.
		If the ``explicit`` and ``listfile`` arguments are None/empty both components of the tuple
		will be None as well.
	"""
	if explicit:
		paths = list(map(Path, explicit))
		paths_str = list(map(str, paths))

	elif listfile is not None:
		lines = list(read_lines(listfile, skip_empty=True))
		paths = [Path(listfile_dir) / line for line in lines]
		paths_str = lines

	else:
		return None, None

	files = SequenceFile.from_paths(paths, 'fasta', 'auto')
	ids = [get_file_id(f, strip_dir, strip_ext) for f in paths_str]

	return ids, files


def warn_duplicate_file_ids(ids: List[str], template: str):
	"""Print a warning message if duplicate file IDs are present.

	Parameters
	----------
	ids
		List of file ID strings, such as from :func:`.get_sequence_files`.
	template
		Message template. May contain formatting placeholders for ``ids`` (comma-delimited string of
		duplicated IDs), ``id`` (first duplicated ID), and ``n`` (number of duplicated IDs).
	"""
	counts = Counter(ids)
	duplicates = [id_ for id_, count in counts.items() if count > 1]
	if duplicates:
		msg = template.format(id=duplicates[0], ids=', '.join(duplicates), n=len(duplicates))
		click.echo(msg, err=True)


################################################################################
# Click introspection
################################################################################

def params_by_name(cmd: click.Command, names: Optional[Iterable[str]]=None):
	"""Get parameters of click command by name.

	Parameters
	----------
	cmd
	names
		Names of specific parameters to get.

	Returns
	-------
	Union[Dict[str, click.Parameter], List[click.Parameter]]
		Parameters with given in ``names`` argument if not None, otherwise a dictionary containing
		all of the command's parameters keyed by name.
	"""
	by_name = {param.name: param for param in cmd.params}
	if names is None:
		return by_name
	else:
		return [by_name[name] for name in names]

def check_params_group(ctx: click.Context, names: Iterable[str], exclusive: bool, required: bool):
	"""Check for the presence of the given parameter values and raise an informative error if needed.

	Parameters
	----------
	ctx
	names
		Parameter names.
	exclusive
		No more than one of the parameters may be present.
	required
		At least one of the parameters must be present.

	Raises
	------
	click.ClickException
	"""
	nfound = sum(bool(ctx.params[name]) for name in names)

	if exclusive and nfound > 1:
		params = params_by_name(ctx.command, names)
		plist = join_list_human(map(param_name_human, params), 'and')
		raise click.ClickException(f'{plist} are mutually exclusive')

	if required and nfound == 0:
		params = params_by_name(ctx.command, names)
		plist = join_list_human(map(param_name_human, params), 'or')
		raise click.ClickException(f'One of {plist} is required')


def param_name_human(param: click.Parameter) -> str:
	"""Get the name/metavar of the given parameter as it appears in the auto-generated help output."""
	if isinstance(param, click.Option):
		# return param.opts[0]
		return '/'.join(param.opts)
	if isinstance(param, click.Argument):
		if param.metavar is not None:
			return param.metavar.rstrip('.')  # Remove ellipsis
		else:
			return param.opts[0].upper()
	raise TypeError(f'Expected click.Parameter, got {type(param)}')


################################################################################
# Misc
################################################################################

def print_table(rows: Sequence[Sequence], colsep: str=' ', left: str='', right: str=''):
	"""Print a basic table."""

	echo = lambda s: click.echo(s, nl=False)

	rows = [list(map(str, row)) for row in rows]
	ncol = max(map(len, rows))

	widths = [0] * ncol
	for row in rows:
		for i, val in enumerate(row):
			widths[i] = max(widths[i], len(val))

	for row in rows:
		echo(left)

		for i, val in enumerate(row):
			echo(val.ljust(widths[i]))

			if i < ncol - 1:
				echo(colsep)

		echo(right)
		echo('\n')
