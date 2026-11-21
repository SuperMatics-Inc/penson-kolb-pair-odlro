# pyright: reportMissingImports=false
#!/usr/bin/env python3
"""
mutants.py — Adversarial mutation test suite for PAIR-LRO.

Verifies that the verification guards in replay.py have real teeth:
  Mutant 1: Replace pi-flux hopping on the 2x2 plaquette with standard 0-flux real hopping
            -> 14x14 Majorana RP Gram matrix positivity MUST FAIL (min_eig < -0.14).
  Mutant 2: Flip the pair-coupling sign J -> -J under the FM unitary (wrong-sign pair hopping
            across the reflection cut) -> 14x14 Majorana RP Gram matrix MUST FAIL (min_eig < -0.15).
  Mutant 3: Use Koma's uncorrected factor-of-2 Gaussian domination formula 1/(2*beta*g*E_{p+Q})
            on the 4-site ring (both at U=0 and U=1) -> Duhamel bound check MUST FAIL (ratio > 1.30 > 1).
  Mutant 4: Raise temperature T0 above the physical XY critical ceiling (e.g. T0/g = 2.5)
            -> Certifier MUST REFUSE to certify (m2_lo < 0).
  Mutant 5: Excessive single-electron hopping |kappa|/g = 1.50 at U=0
            -> Certifier MUST REFUSE to certify (m2_lo < 0).

Exits 0 iff every mutant is caught as expected.
"""

import math
import os
import sys
import numpy as np
import mpmath as mp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import replay  # noqa: E402

print("\n" + "=" * 72)
print("RUNNING ADVERSARIAL MUTATION SUITE (mutants.py)")
print("=" * 72)

# Mutant 1: 0-flux real hopping must violate Majorana RP on 2x2 plaquette
_, eig_0flux_fm = replay.test_rp_plaquette(
    replay.H_plaq_pair_FM + replay.H_plaq_U + replay.H_plaq_real0,
    replay.U_ph_odd @ replay.U_phase_odd,
    beta_val=2.0,
)
assert eig_0flux_fm < -0.1, f"Mutant 1 survived! eig_0flux_fm={eig_0flux_fm}"
print(f"  [CAUGHT] Mutant 1 (0-flux real hopping breaks Majorana RP): min_eig = {eig_0flux_fm:.6f} < 0")

# Mutant 2: Flipping pair-coupling sign J -> -J under U_ph * U_phase must violate Majorana RP
_, eig_wrong_J = replay.test_rp_plaquette(
    -replay.H_plaq_pair_FM + replay.H_plaq_U + replay.H_plaq_pi_real,
    replay.U_ph_odd @ replay.U_phase_odd,
    beta_val=2.0,
)
assert eig_wrong_J < -0.1, f"Mutant 2 survived! eig_wrong_J={eig_wrong_J}"
print(f"  [CAUGHT] Mutant 2 (Flipping J -> -J breaks Majorana RP):    min_eig = {eig_wrong_J:.6f} < 0")

# Mutant 3: Koma's 1/(2*beta*g*E_{p+Q}) formula must be violated by exact Duhamel b_p on 4-site ring
for label_m3, ev_m3, evec_m3, Q_m3 in [
    ("Koma U=0 AF (kap=0.0)", replay.evals4_koma0, replay.evecs4_koma0, math.pi),
    ("Koma U=0 AF (kap=0.1)", replay.evals4_koma, replay.evecs4_koma, math.pi),
    ("FM + U=1    (kap=0.1)", replay.evals4, replay.evecs4, 0.0),
]:
    max_koma_ratio_m3 = 0.0
    for beta_val in [1.0, 2.0]:
        w = np.exp(-beta_val * ev_m3)
        Z = float(np.sum(w))
        dE = ev_m3[None, :] - ev_m3[:, None]
        dw = w[:, None] - w[None, :]
        K = np.where(np.abs(dE) > 1e-10, dw / np.where(np.abs(dE) > 1e-10, beta_val * dE, 1.0), w[:, None])
        for k_idx in range(4):
            p = 2.0 * math.pi * k_idx / 4.0
            EpQ = 1.0 - math.cos(p + Q_m3)
            if EpQ < 1e-9:
                continue
            Gp = replay.Gp_list[k_idx]
            Ad = evec_m3.conj().T @ Gp.conj().T @ evec_m3
            Be = evec_m3.conj().T @ Gp @ evec_m3
            b_val = float(np.sum(K * Ad * Be.T).real) / Z
            koma_uncorrected_bound = 1.0 / (2.0 * beta_val * 1.0 * EpQ)
            max_koma_ratio_m3 = max(max_koma_ratio_m3, b_val / koma_uncorrected_bound)
    assert max_koma_ratio_m3 > 1.30, f"Mutant 3 ({label_m3}) survived! max_koma_ratio_m3={max_koma_ratio_m3}"
    print(
        f"  [CAUGHT] Mutant 3 ({label_m3} Koma 1/(2*beta*g*E) slip): max b_p / [1/(2*beta*g*E)] = {max_koma_ratio_m3:.6f} > 1.0"
    )

# Mutant 4: Excessive temperature T0/g = 2.5 at U/g = 10 must NOT be certified
m2_mutant4 = replay.certify_all_T_le_T0(2500, 1000, 10, 1, 0, 1, 1, 1, replay.Jsq_1)
assert m2_mutant4 < 0, f"Mutant 4 survived! m2_mutant4={m2_mutant4}"
print(
    f"  [CAUGHT] Mutant 4 (Unphysical high T0/g = 2.50): certifier returned m_0^2 >= {mp.nstr(m2_mutant4, 5)} < 0 (REJECTED)"
)

# Mutant 5: Excessive hopping |kappa|/g = 1.5 at U/g = 0, T0/g = 0.16 must NOT be certified
m2_mutant5 = replay.certify_all_T_le_T0(160, 1000, 0, 1, 150, 100, 1, 2, replay.Jsq_half)
assert m2_mutant5 < 0, f"Mutant 5 survived! m2_mutant5={m2_mutant5}"
print(
    f"  [CAUGHT] Mutant 5 (Excessive hopping |kappa|/g = 1.50 at U=0): certifier returned m_0^2 >= {mp.nstr(m2_mutant5, 5)} < 0 (REJECTED)"
)

print("\nALL 5 ADVERSARIAL MUTANTS CAUGHT SUCCESSFULLY.")
