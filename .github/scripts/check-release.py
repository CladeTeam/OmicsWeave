"""Keep Python and Rust package versions consistent with an optional release tag."""

import argparse
import re
import tomllib
from pathlib import Path

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--tag")
args = parser.parse_args()
root = Path(__file__).resolve().parents[2]
python = tomllib.loads((root / "pyproject.toml").read_text())["project"]["version"]
rust = tomllib.loads((root / "Cargo.toml").read_text())["workspace"]["package"]["version"]
if python != rust:
    parser.error(f"Python version {python!r} differs from Rust version {rust!r}")
if args.tag is not None:
    if not re.fullmatch(r"v\d+\.\d+\.\d+", args.tag):
        parser.error("A stable release tag must have the form vMAJOR.MINOR.PATCH")
    if args.tag != f"v{python}":
        parser.error(f"Tag {args.tag!r} does not match package version v{python}")
print(f"Package versions agree: {python}")
