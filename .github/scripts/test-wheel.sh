#!/usr/bin/env bash
set -euo pipefail

WHEEL_DIR=$(cd "$1" && pwd)
shift
PROJECT_DIR=$(cd "$(dirname "$0")/../.." && pwd)
TEST_DIR=$(mktemp -d "${RUNNER_TEMP:-${TMPDIR:-/tmp}}/omicsweave-wheel.XXXXXX")
trap 'rm -rf "$TEST_DIR"' EXIT
cp -R "$PROJECT_DIR/src/tests" "$TEST_DIR/tests"
cd "$TEST_DIR"

for version in "$@"; do
    uv python install "$version"
    uv venv --managed-python --python "$version" "env-$version"
    PYTHON="$TEST_DIR/env-$version/bin/python"
    uv pip install --python "$PYTHON" "$WHEEL_DIR"/*.whl 'pytest>=8' "${NUMPY_SPEC:-numpy}"
    uv pip check --python "$PYTHON"
    "$PYTHON" -I -c 'import importlib.metadata as m; import numpy as np; import omicsweave as ow; assert ow.__version__ == m.version("omicsweave"); print(ow.__file__, ow.__version__, "NumPy", np.__version__)'
    "$TEST_DIR/env-$version/bin/ow" --version
    "$PYTHON" -m pytest -ra tests
done
