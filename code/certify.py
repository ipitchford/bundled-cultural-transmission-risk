"""Small helpers shared by the certificate scripts.

Checks are explicit (never `assert`), so they also run under `python -O`.
A script records every check, writes its status string from the actual
results and ends with exit status 1 if any check failed.
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_RESULTS = ROOT / "results"
DEFAULT_FIGURES = ROOT / "figures"


class Checker:
    """Record named pass/fail checks without relying on `assert`."""

    def __init__(self) -> None:
        self.count = 0
        self.failures: list[dict] = []

    def check(self, condition, label: str, **detail) -> bool:
        self.count += 1
        try:
            ok = bool(condition)
        except Exception:  # e.g. ambiguous array truth value
            ok = False
        if not ok:
            record = {"label": label}
            record.update({k: _jsonable(v) for k, v in detail.items()})
            self.failures.append(record)
            print(f"CHECK FAILED: {label} {record}", file=sys.stderr, flush=True)
        return ok

    def close(self, label: str, got, expected, tol: float, relative: bool = False, **detail) -> bool:
        got_f, exp_f = float(got), float(expected)
        err = abs(got_f - exp_f)
        if relative:
            err /= max(abs(exp_f), 1e-300)
        return self.check(math.isfinite(err) and err <= tol, label, got=got_f,
                          expected=exp_f, error=err, tol=tol, **detail)

    @property
    def ok(self) -> bool:
        return not self.failures

    def status(self) -> str:
        if self.ok:
            return f"all {self.count} checks passed"
        return f"FAILED: {len(self.failures)} of {self.count} checks failed"

    def summary(self) -> dict:
        return {"status": self.status(), "checks_run": self.count,
                "checks_failed": len(self.failures), "failures": self.failures}

    def finish(self) -> None:
        print(self.status(), flush=True)
        if self.failures:
            sys.exit(1)


def _jsonable(x):
    try:
        json.dumps(x)
        return x
    except TypeError:
        return str(x)


def output_dir(arg: str | None, default: Path) -> Path:
    """Resolve --out. Writing into the default directory may refresh files."""
    path = default if arg is None else Path(arg).expanduser().resolve()
    path.mkdir(parents=True, exist_ok=True)
    return path


def safe_write(path: Path, data, default: Path, binary: bool = False) -> Path:
    """Write a result file; refuse to overwrite outside the default directory."""
    path = Path(path)
    if path.exists() and path.resolve().parent != Path(default).resolve():
        raise SystemExit(f"Refusing to overwrite existing file outside the default "
                         f"results directory: {path}")
    if binary:
        path.write_bytes(data)
    else:
        path.write_text(data)
    return path


def check_writable(path: Path, default: Path) -> Path:
    """Fail early (before long computations) if a later write would be refused."""
    path = Path(path)
    if path.exists() and path.resolve().parent != Path(default).resolve():
        raise SystemExit(f"Refusing to overwrite existing file outside the default "
                         f"results directory: {path}")
    return path
