# Explicit positive-temperature pair long-range order in Penson-Kolb-Hubbard models
Research record: [Explicit positive-temperature pair long-range order in Penson-Kolb-Hubbard models](https://supermatics.io/research/penson-kolb-pair-odlro).

Davis Rens, SuperMatics Inc., [ORCID 0000-0001-5839-5197](https://orcid.org/0000-0001-5839-5197)

| | |
| --- | --- |
| Preprint | [10.5281/zenodo.24428099](https://doi.org/10.5281/zenodo.24428099) |
| Replay supplement | [10.5281/zenodo.24428101](https://doi.org/10.5281/zenodo.24428101) |
| Version | 1.0.0, 21 November 2026 |
| Manuscript | [`paper/main.pdf`](paper/main.pdf), 24 pages |

This repository is the tagged release tree of that preprint. It holds the reviewed manuscript
source, the released PDF and the complete replay supplement, pinned byte for byte through
[`MANIFEST.sha256`](MANIFEST.sha256) and [`replay/SHA256SUMS`](replay/SHA256SUMS).

## What the paper does

The paper proves explicit, constructive lower bounds on onsite pair off-diagonal long-range order
for three-dimensional Penson-Kolb-Hubbard models on the full fermion Fock space, at every
temperature `0 < k_B T <= T_0`, together with a uniform condensate floor
`inf_{0 < T <= T_0} m_0^2(T) > 0` and explicit floors `T_0` in the thermodynamic limit. It covers
two regimes and says where the method stops.

1. **Zero single-electron hopping.** With `t = 0`, pair-hopping amplitude `J > 0` and on-site
   coupling `U >= -0.5J`, so that weak on-site repulsion is covered as well as attraction, the
   certified floors run from `T_0/J = 0.155` at `U = -0.5J` to `T_0/J = 1.170` at `U = 100J`. At
   `t = 0` there is no single-electron hopping, so the flux plays no role.
2. **Pi-flux single-electron hopping.** With `|kappa|/g <= 0.50` and on-site attraction `U >= 0`,
   Majorana reflection positivity together with a Kubo-Kishi weighted infrared sum rule gives the
   same bound, with certified floors `T_0/g` up to `0.590` at `|\kappa|/g = 0.10`, `U/g = 10`.
3. **One certification for both regimes.** Both are reached by a single finite-temperature
   variational and entropic argument. It combines Majorana reflection positivity with arbitrary
   on-site Hubbard interactions, a corrected Gaussian domination infrared bound, a
   vacancy-controlled Cauchy-Schwarz bound on single-electron hopping, a KLS-weighted infrared sum
   rule, a vacancy-conditioned von Neumann entropy bound, an exact four-state dimer trial state and
   an exact seven-site star operator-norm bound.
4. **A factor of two in a printed constant.** The Gaussian domination constant used in the cited
   work is off by a factor of two. An exact four-site ring comparison against the printed inequality
   gives ratio `1.329280 > 1`, and the `2 x 2` plaquette Majorana Gram matrix has
   `lambda_min = -0.146444` for zero-flux hopping.
5. **Finite arithmetic, verified by outward-rounded interval arithmetic.** Five families of
   Brillouin-zone integrals are enclosed at 50 digits with outward-rounded interval arithmetic, the
   fifteen rows of the two floor tables are certified the same way, and the finite operator,
   double-commutator, dimer, star, ring and plaquette claims are checked by exact
   finite-Fock-space algebra and exact diagonalisation. Five semantic mutants, one per claim
   class, are each rejected.

## What the paper does not claim

Positive-temperature pair long-range order for zero-flux real single-electron hopping at `t != 0`
remains open, and the paper documents three routes that do not reach it: Majorana reflection
positivity fails there, spin reflection positivity does not reach spatial pair infrared bounds, and
norm or contour expansions about `t = 0` do not converge. No material is claimed and no
superconductivity is claimed. No phase transition and no upper bound on a critical temperature is
proved or excluded, and no statement is made about a measured transition temperature anywhere in
the paper. The interval enclosures and the exact diagonalisations check finite arithmetic; they do
not formalise the analytic or thermodynamic arguments that carry the theorems.

## Contents

| Path | What it is |
| --- | --- |
| `paper/main.pdf` | The released manuscript, 24 pages. |
| `paper/main.tex` | The LaTeX source that produces that PDF. |
| `paper/refs.bib` | The bibliography. |
| `BUILD.md` | How the PDF is built, with the compiler version, the fixed epoch and the observed log counts. |
| `replay/README.md` | What each replay script checks, how to run it and its byte custody. |
| `replay/replay.py` | The Brillouin-zone interval enclosures, the fifteen floor rows, the exact finite-Fock-space identities, the exact diagonalisations, the ring and plaquette comparisons. |
| `replay/mutants.py` | Five semantic mutations, one per claim class, each of which must be rejected. |
| `replay/run_all.py` | Standard-library driver that runs both scripts and reports each exit code. |
| `replay/SHA256SUMS` | Pins for the four other files of the replay directory. |
| `MANIFEST.sha256` | Pins for every other tracked file of this repository. |
| `CITATION.cff` | Citation metadata for the preprint and the supplement. |
| `LICENSE` | Apache-2.0 for the code, CC-BY-4.0 for the manuscript, as stated in its first lines. |
| `.github/workflows/ci.yml` | The two hosted jobs described below. |

`paper/main.tex` differs from the admitted reviewed source only on its two DOI macros,
`\paperdoi` and `\suppdoi`, which carry the archival identifiers of this release. `paper/refs.bib`
is unchanged. `replay/` is unchanged, and its own `SHA256SUMS` still verifies.

## Reproduce

### The replay

Python 3.12 or later. The replay scripts import only NumPy and mpmath, at the verified versions:

```sh
python -m pip install numpy==2.5.5 mpmath==1.4.1
```

From the `replay/` directory of an unpacked release:

```sh
cd replay
python -B replay.py
python -B mutants.py
python -B run_all.py
```

Every command exits `0`. The closing lines are `ALL REPLAY CHECKS PASSED` and
`ALL 5 ADVERSARIAL MUTANTS CAUGHT SUCCESSFULLY`, and `run_all.py` prints
`SUMMARY: 2/2 scripts passed in ...s` and exits nonzero if either script fails. Timings vary with
the machine and its load; the closing text does not. `run_all.py` takes about a minute on a single
laptop core. [`replay/README.md`](replay/README.md) maps every check to the statement it backs and
lists the observed final lines.

### The paper

```sh
cd paper
SOURCE_DATE_EPOCH=1795132800 tectonic -X compile --keep-logs main.tex
```

Tectonic 0.17.0 with that fixed epoch reproduces the released PDF byte for byte (SHA-256
`35b60cbbe38801a5a2e38f4dee8cf9b1ec025a7854bd11e1df770cf107613b57`); the log has no overfull box,
no undefined reference, no warning and three underfull boxes documented in
[`paper/BUILD.md`](paper/BUILD.md).

### Continuous integration

[`.github/workflows/ci.yml`](.github/workflows/ci.yml) has two jobs. The `replay` job verifies
`MANIFEST.sha256` and `replay/SHA256SUMS`, installs the pinned NumPy, SciPy, mpmath and sympy,
and runs `run_all.py`, that is `replay.py` followed by `mutants.py`. The `paper` job builds the
manuscript with the epoch above, requires zero overfull boxes, zero undefined references, zero
warnings and exactly 24 pages, and requires the paper and supplement DOIs in the built PDF.

## Byte custody

Every released file is a copy of reviewed bytes, and `MANIFEST.sha256` plus `replay/SHA256SUMS`
pin all of them except the two manifests themselves:

```sh
sha256sum -c MANIFEST.sha256
(cd replay && sha256sum -c SHA256SUMS)
```

The deposited supplement archive is `git archive` of the `v1.0.0` tag, so its members are the
files listed by `git ls-files` at that commit, unchanged.

## Licences

The manuscript (`paper/main.tex`, `paper/refs.bib`, `paper/main.pdf`) is licensed under Creative
Commons Attribution 4.0 International (CC-BY-4.0). The replay scripts and all remaining files are
licensed under the Apache License 2.0. The first lines of [`LICENSE`](LICENSE) state that split;
the full Apache-2.0 text follows.

## Review and AI assistance

AI agents generated proofs, certificates, code, and manuscript drafts and edits under D. Rens'
direction and review. A separate AI review of the repaired draft returned PASS with seven fixes
applied. This preprint has not been peer reviewed.
