# Developer Guide

## Releases and publishing

Distributions are built and published to PyPI by the `Release` GitHub Actions workflow
(`.github/workflows/release.yml`).


### Workflow overview

The workflow has two triggers:

- **Publishing a GitHub release** builds everything, attaches the wheels and sdist to the release,
  and publishes them to PyPI. Note that this includes releases marked as prereleases.
- **Manual dispatch** (Actions tab → Release → "Run workflow") builds everything and uploads the
  files as workflow artifacts. The `publish` input controls where they go next: `none` (the
  default) stops there, and `testpypi` also publishes to TestPyPI.

Jobs:

| Job                 | Runs on                     | Does                                                        |
|---------------------|-----------------------------|-------------------------------------------------------------|
| `build-wheels`      | always                      | Builds wheels with cibuildwheel, uploads `wheels` artifact  |
| `build-sdist`       | always                      | Runs `build-sdist` and `check-dist` Pixi tasks, uploads `sdist` artifact |
| `attach-to-release` | release                     | Attaches all wheels and the sdist to the GitHub release     |
| `publish-pypi`      | release                     | Publishes to PyPI                                           |
| `publish-testpypi`  | dispatch with `testpypi`    | Publishes to TestPyPI                                       |


### Authentication with PyPI (trusted publishing)

The publish jobs use [trusted publishing](https://docs.pypi.org/trusted-publishers/): GitHub
issues a short-lived OIDC token to the job, and PyPI exchanges it for an upload token. No API
tokens or repository secrets are needed. This requires three things.

1. **GitHub environments.** In the repository settings (Settings → Environments), create
   environments named `pypi` and `testpypi`. Adding required reviewers to `pypi` makes every PyPI
   upload wait for manual approval.

2. **Trusted publishers on PyPI and TestPyPI.** In the `gambit` project settings on
   [pypi.org](https://pypi.org/manage/project/gambit/settings/publishing/) and
   [test.pypi.org](https://test.pypi.org/manage/project/gambit/settings/publishing/), add a GitHub
   publisher. If the project doesn't exist yet on an index, add a "pending publisher" from your
   account's publishing page instead.

   | Field             | pypi.org      | test.pypi.org |
   |-------------------|---------------|---------------|
   | Owner             | `jlumpe`      | `jlumpe`      |
   | Repository        | `gambit`      | `gambit`      |
   | Workflow filename | `release.yml` | `release.yml` |
   | Environment       | `pypi`        | `testpypi`    |

   The workflow filename and environment name must match the workflow exactly. Renaming either
   breaks publishing until the publisher entry is updated.

3. **Job permissions.** The publish jobs need `id-token: write` to request the OIDC token. This is
   already set in the workflow, and the top-level permissions are kept to `contents: read`.

TestPyPI is a separate service from PyPI with its own accounts, so the TestPyPI publisher has to be
set up from a test.pypi.org account.


### Making a release

1. Set `version` in `[project]` in `pyproject.toml` and add an entry to `CHANGELOG.md`.
2. Optionally check the packaging locally:

   ```bash
   pixi run -e dist build-sdist
   pixi run -e dist check-dist
   ```

3. Optionally do a TestPyPI dry run (see below).
4. Merge to `main`, tag the release (e.g. `v1.2.0`), and publish a GitHub release from the tag.
   This triggers the workflow, which attaches the files to the release and publishes to PyPI.


### TestPyPI dry run

1. Set a dev version (e.g. `1.2.0.dev1`) on a branch and push it.
2. Run the workflow manually on that branch with `publish` set to `testpypi`.
3. Install from TestPyPI, pulling dependencies from the real PyPI:

   ```bash
   pip install -i https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple/ gambit==1.2.0.dev1
   gambit --help
   ```

TestPyPI (like PyPI) never accepts the same version twice, even after deleting it, so each retry
needs a new version (`dev2`, `dev3`, ...).
