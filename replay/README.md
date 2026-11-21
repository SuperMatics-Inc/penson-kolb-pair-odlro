# Replay supplement

Accompanying paper: *Explicit positive-temperature pair long-range order in Penson-Kolb-Hubbard models*

Author: Davis Rens, SuperMatics Inc., ORCID 0000-0001-5839-5197

Two standalone Python scripts (`replay.py` and `mutants.py`) and a standard-library
driver (`run_all.py`) reproduce every finite-arithmetic, outward-rounded interval
enclosure, exact finite-Fock-space algebraic identity, small-cluster exact
diagonalization check, and adversarial mutation check referenced in the
accompanying paper. Both `replay.py` and `mutants.py` are byte-for-byte copies
of the reviewed scripts (`mutants.py` imports `replay` from this directory, and
neither script reads any file outside this directory). `SHA256SUMS` in this
directory lists the SHA-256 digest of every other file in this directory, and the
digests of the three Python files are repeated at the end of this file.

| File | Role | Primary checks |
|---|---|---|
| `replay.py` | Main verification and interval-bound replay | 50-digit outward interval enclosures for 3D Brillouin-zone integrals and all 15 `T_0` rows in Tables 1 and 2, exact finite-Fock-space algebraic and double-commutator identities, 4-state dimer ED, 128-state 7-site star norm, 4-site ring Duhamel factor-of-2 check, 2x2 plaquette 14x14 Majorana RP Gram matrix test, and literal formula/table check |
| `mutants.py` | Adversarial mutation suite | 5 semantic mutations (0-flux real hopping on 2x2 plaquette, flipped pair coupling `J -> -J`, uncorrected Koma factor-of-2 Gaussian domination on 4-site ring, `T_0/g = 2.50` above XY ceiling, and `|kappa|/g = 1.50` at `U = 0`) |
| `run_all.py` | Combined driver | Runs `replay.py` then `mutants.py` via `sys.executable -B`, prints each exit code, and exits nonzero if either fails |

## Script descriptions

`replay.py` is a standalone Python verification script (`numpy`, `mpmath`) that
reproduces all 50-digit interval enclosures with strict outward directed
rounding (`I_W`, `I_3`, `W_3(alpha)`, `J_sq(alpha)`, `I_sin`), all 15 explicit
`T_0` rows in Tables 1 and 2 with strict directed rounding down, all exact
finite-Fock-space operator and on-site/cross-site double-commutator identities,
the 4-state dimer exact-diagonalization verification of `F_dimer_pair_upper`,
the 128-state 7-site star norm verification (`2*sqrt(3)`), the 4-site ring
Duhamel factor-of-2 check on both Koma's `U = 0` antiferromagnetic model and
the `U = 1` ferromagnetic model, the 2x2 plaquette 14x14 Majorana
reflection-positivity Gram matrix test (confirming reflection positivity holds
for pi-flux hopping and fails for 0-flux real hopping with minimum eigenvalue
`-0.146444 < 0`), and the Section 6 literal formula and table self-check
against the proof's Table 1 and Table 2 bounds.

`mutants.py` is an adversarial mutation test suite confirming that all 5
semantic mutations are caught and rejected:
(1) replacing pi-flux hopping with 0-flux real hopping on the 2x2 plaquette
(breaking Majorana reflection positivity, `min_eig < -0.14`);
(2) flipping ferromagnetic `J -> -J` on the 2x2 plaquette (breaking Majorana
reflection positivity, `min_eig < -0.15`);
(3) restoring Koma's uncorrected factor of 2 in Gaussian domination on the
4-site ring for both Koma's `U = 0` antiferromagnetic model (`kappa = 0, 0.1`)
and the `U = 1` ferromagnetic model (`kappa = 0.1`), where the exact Duhamel
bound exceeds the uncorrected formula with ratio `> 1.30 > 1`;
(4) raising `T_0/g = 2.50` at `U/g = 10` above the XY critical ceiling
(rejected by the interval bound evaluator with negative lower bound); and
(5) increasing `|kappa|/g = 1.50` at `U = 0`, `T_0/g = 0.16` beyond the
positive-bound window (rejected with negative lower bound).

`run_all.py` is a small standard-library-only driver that executes `replay.py`
and then `mutants.py` via `sys.executable -B`. It prints a `PASS` or `FAIL`
line with each script's exit code (`rc=...`), wall time, and final output line,
and exits nonzero if either script fails.

## Requirements and commands

The scripts require Python 3.12 or newer with `numpy` and `mpmath`. They were
verified with Python 3.13.9, `numpy` 2.5.5, and `mpmath` 1.4.1. From this
directory (`supp/`) in a POSIX shell with `uv`, run the full driver with:

```sh
uv run --no-project --with numpy --with mpmath python run_all.py
```

To run each script individually:

```sh
uv run --no-project --with numpy --with mpmath python replay.py
uv run --no-project --with numpy --with mpmath python mutants.py
```

In an existing Python 3.12+ environment with `numpy` and `mpmath` installed,
replace the `uv run ...` prefix with `python`. Setting `OPENBLAS_NUM_THREADS=1`,
`OMP_NUM_THREADS=1`, and `MKL_NUM_THREADS=1` avoids BLAS thread contention, and
`PYTHONDONTWRITEBYTECODE=1` prevents Python from creating `__pycache__`
artifacts (`run_all.py` sets these defaults automatically for child processes).

## Expected output and approximate run times

Both scripts and `run_all.py` exit with status 0. Observed final output lines
and laptop wall times from `run_all.py` are:

- `replay.py` (~24.7 s, `rc=0`): `ALL REPLAY CHECKS PASSED in 24.44s.`
- `mutants.py` (~27.8 s, `rc=0`): `ALL 5 ADVERSARIAL MUTANTS CAUGHT SUCCESSFULLY.`
- `run_all.py` (~52.5 s total, `rc=0`): `SUMMARY: 2/2 scripts passed in 52.50s`

## What these finite checks do not establish

These finite checks do not by themselves prove the thermodynamic-limit theorem
in the paper. Passing the 50-digit outward interval evaluations, finite-Fock-space
algebraic identities, small-cluster exact diagonalizations, and mutation tests
verifies finite arithmetic and finite-volume witnesses only; the thermodynamic
lower bounds follow from the analytical arguments in the paper. All bounds
concern mathematical pair correlations of idealized lattice fermion Hamiltonians,
and positive-temperature pair off-diagonal long-range order for 0-flux real
single-electron hopping (`t != 0`) remains open.

## Licensing and AI disclosure

Code in this supplement (`replay.py`, `mutants.py`, `run_all.py`) is licensed
under the Apache License, Version 2.0 (Apache-2.0). The accompanying paper is
distributed separately under CC-BY-4.0.

AI agents generated proofs, certificates, code and manuscript drafts and edits
under the author's direction and review.

## SHA-256 of the Python scripts in this directory

```
51fb20bfd974f5736bd730c02f561ac46469f8f52c4b855b779894601883387c  replay.py
db3ffad4024cff3a8fcb51dcb1e5b5a12fbc8aac8269879ab8a3cbb67056c03d  mutants.py
407ed369f4f2770d0c6972db66353decc138e47abeb4ef5d4a574e35e5150644  run_all.py
```

`SHA256SUMS` covers every file in `supp/` except `SHA256SUMS` itself (including
this `README.md`) and can be verified from `supp/` with `sha256sum -c SHA256SUMS`.
