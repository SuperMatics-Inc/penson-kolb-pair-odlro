#!/usr/bin/env python3
"""Standard-library-only driver for the Penson-Kolb-Hubbard ODLRO reproducibility supplement."""

from __future__ import annotations

import os
from pathlib import Path
import subprocess
import sys
import time


SCRIPTS: list[str] = [
    "replay.py",
    "mutants.py",
]


def main() -> int:
    supp_dir = Path(__file__).resolve().parent
    env = dict(os.environ)
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    env.setdefault("OPENBLAS_NUM_THREADS", "1")
    env.setdefault("OMP_NUM_THREADS", "1")
    env.setdefault("MKL_NUM_THREADS", "1")

    t_total_start = time.monotonic()
    all_ok = True
    results: list[tuple[str, int, float]] = []

    for script in SCRIPTS:
        cmd = [sys.executable, "-B", str(supp_dir / script)]
        t0 = time.monotonic()
        proc = subprocess.run(
            cmd,
            cwd=supp_dir,
            env=env,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            encoding="utf-8",
            errors="replace",
        )
        dt = time.monotonic() - t0
        lines = [ln for ln in proc.stdout.splitlines() if ln.strip()]
        tail = lines[-1] if lines else "(no output)"
        if proc.returncode == 0:
            print(f"PASS  {script}  rc=0  {dt:.2f}s  | {tail}", flush=True)
        else:
            all_ok = False
            print(f"FAIL  {script}  rc={proc.returncode}  {dt:.2f}s  | {tail}", flush=True)
            for extra in lines[-10:]:
                print(f"      {extra}", flush=True)
        results.append((script, proc.returncode, dt))

    total_dt = time.monotonic() - t_total_start
    passed = sum(1 for _, rc, _ in results if rc == 0)
    print(f"SUMMARY: {passed}/{len(results)} scripts passed in {total_dt:.2f}s", flush=True)
    return 0 if all_ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
