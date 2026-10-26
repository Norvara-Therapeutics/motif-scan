"""Build a wheel and test the *installed* package, not the source tree.

Running the suite from a checkout proves the source works. It does not prove
the wheel does: a module left out of the package, a missing entry point or a
data file that never made it into the distribution all pass a source-tree run
and fail for whoever installs it.

This builds the distributions, installs the wheel into a throwaway virtual
environment with ``PYTHONPATH`` stripped so nothing resolves back to ``src/``,
and runs the full suite plus a command-line smoke test from a directory that is
not the repository.

    python scripts/check_wheel.py
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
import tempfile
import venv
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def run(command: list[str], cwd: Path, env: dict[str, str]) -> subprocess.CompletedProcess:
    print(f"$ {' '.join(str(part) for part in command)}")
    return subprocess.run(command, cwd=cwd, env=env, text=True)


def main() -> int:
    clean = {key: value for key, value in os.environ.items() if key != "PYTHONPATH"}

    dist = ROOT / "dist"
    if dist.exists():
        shutil.rmtree(dist)
    built = run(
        [sys.executable, "-m", "build", "--no-isolation"], ROOT, clean
    )
    if built.returncode:
        print("build failed", file=sys.stderr)
        return 1

    wheels = sorted(dist.glob("*.whl"))
    if len(wheels) != 1:
        print(f"expected exactly one wheel, found {len(wheels)}", file=sys.stderr)
        return 1
    wheel = wheels[0]

    with tempfile.TemporaryDirectory(prefix="motif-scan-wheel-") as holding:
        home = Path(holding)
        environment = home / "venv"
        venv.create(environment, with_pip=True)
        binaries = environment / ("Scripts" if os.name == "nt" else "bin")
        python = binaries / "python"

        install = run(
            [str(python), "-m", "pip", "install", "--quiet", f"{wheel}[dev]"],
            home,
            clean,
        )
        if install.returncode:
            print("install failed", file=sys.stderr)
            return 1

        # Tests and examples are not installed by the wheel, so they are run
        # from the checkout while the package itself resolves to the venv.
        suite = run(
            [str(python), "-m", "pytest", "-q", str(ROOT / "tests")],
            ROOT,
            clean,
        )
        if suite.returncode:
            print("installed package fails its own test suite", file=sys.stderr)
            return 1

        checks: list[list[str]] = [
            [str(binaries / "motif-scan"), "--version"],
            [str(binaries / "motif-scan"), "motifs", "--type", "protein"],
            [str(binaries / "motif-scan"), "codes", "--type", "nucleotide"],
            [
                str(binaries / "motif-scan"),
                "scan",
                "-i",
                str(ROOT / "examples" / "aav2.fasta"),
                "-m",
                "TATA,polyA,AAV-RBE",
                "--mismatches",
                "1",
                "-f",
                "json",
            ],
            [
                str(binaries / "motif-scan"),
                "scan",
                "-i",
                str(ROOT / "examples" / "aav2_vp1.fasta"),
                "--all",
                "-f",
                "none",
            ],
            [str(python), "-m", "motif_scan", "--version"],
        ]
        for command in checks:
            smoke = run(command, home, clean)
            if smoke.returncode:
                print(f"smoke test failed: {' '.join(command)}", file=sys.stderr)
                return 1

    print("\nwheel verified: installed package passes the suite and the CLI runs")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
