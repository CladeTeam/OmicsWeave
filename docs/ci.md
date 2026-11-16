# CI and PyPI releases

The `CI` workflow runs on pull requests, pushes to `main`, and manual dispatch.
The `CI passed` job is the aggregate status for branch protection.

- Ruff, rustfmt and actionlint check Python, Rust and workflow files.
- Cargo checks all workspace targets on the declared minimum Rust version,
  then runs Clippy and the Rust tests on stable, including the vendored
  reference-tool comparisons. Tests marked `ignore` need external datasets
  or are manual benchmarks and are not part of CI.
- Release wheels are built for Linux x86_64 / AArch64 (manylinux_2_28) and
  macOS Intel / Apple Silicon (macOS 11+). The native engine uses Unix APIs;
  Windows wheels are not built.
- Each wheel is installed into isolated environments outside the checkout
  and tested on CPython 3.12, 3.13, 3.14 and 3.15. The extension uses the
  CPython 3.12 stable ABI, so one wheel serves these Python versions.
  Platform-specific tests skip where appropriate; optional `pysam`
  comparisons skip when it is not installed.
- The source distribution is rebuilt into a wheel and tested separately.
  Wheel and source metadata are checked with Twine. Successful builds are
  retained as Actions artifacts for 14 days.

## Publish a stable release

1. Set the same version in `pyproject.toml` and `[workspace.package]` in
   `Cargo.toml`, refresh `uv.lock` and `Cargo.lock`, and merge the change.
2. Wait for `CI passed` on that commit.
3. Create a tag `vMAJOR.MINOR.PATCH` on that commit and publish its GitHub
   Release. For example, version `0.1.0` requires tag `v0.1.0`.
4. `Publish to PyPI` checks the tag against both package versions, reruns
   the entire CI workflow, and uploads the four tested wheels and source
   distribution after all checks pass.

The upload job uses the `pypi` GitHub environment and its `PYPI_API_TOKEN`
secret. Tokens must never be committed. Upload runs only for a published,
non-prerelease GitHub Release in `CladeTeam/OmicsWeave`. Pushes, tags alone,
pull requests, draft releases, and manual CI dispatch do not publish.
There is deliberately no manual publish dispatch or automatic version bump.

To validate the pipeline without publishing:

```sh
gh workflow run ci.yml --ref main
gh run list --workflow ci.yml
gh run view <run-id> --log-failed
gh run download <run-id> --pattern 'wheel-*' --dir downloaded-wheels
```

The publishing job downloads artifacts from its own run; it does not build
code or receive artifacts from pull-request runs. With API-token publishing,
PyPI attestations are disabled. Re-running an already uploaded release does
not silently skip existing files; inspect any partial upload before retrying.
