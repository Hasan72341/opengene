from __future__ import annotations

import os
from pathlib import Path
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[1]
PACKAGE_DIRS = (
    "core",
    "io",
    "restriction",
    "primers",
    "annotations",
    "alignment",
    "sequencing",
    "cloning",
    "pcr",
    "api",
)


def main() -> int:
    package_paths = [str(ROOT / "packages" / name) for name in PACKAGE_DIRS]
    environment = {**os.environ, "PYTHONPATH": os.pathsep.join(package_paths + [os.environ.get("PYTHONPATH", "")])}
    test_dirs = [ROOT / "packages" / name / "tests" for name in PACKAGE_DIRS]
    test_dirs.append(ROOT / "tests")
    for test_dir in test_dirs:
        sources = "\n".join(path.read_text(encoding="utf-8") for path in test_dir.glob("test_*.py"))
        if "unittest.TestCase" in sources:
            result = subprocess.run(
                [sys.executable, "-m", "unittest", "discover", "-s", str(test_dir), "-v"],
                cwd=ROOT,
                env=environment,
                check=False,
            )
            if result.returncode:
                return result.returncode
        function_runner = """
import pathlib
import runpy
import sys

test_dir = pathlib.Path(sys.argv[1])
failures = []
for path in sorted(test_dir.glob('test_*.py')):
    namespace = runpy.run_path(str(path))
    for name, value in sorted(namespace.items()):
        if name.startswith('test_') and callable(value):
            try:
                value()
            except Exception as error:
                failures.append(f'{path.name}::{name}: {error!r}')
if failures:
    print('\\n'.join(failures), file=sys.stderr)
    raise SystemExit(1)
"""
        function_result = subprocess.run(
            [sys.executable, "-c", function_runner, str(test_dir)],
            cwd=ROOT,
            env=environment,
            check=False,
        )
        if function_result.returncode:
            return function_result.returncode
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
