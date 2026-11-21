# pyright: reportOperatorIssue=false, reportArgumentType=false, reportAttributeAccessIssue=false, reportInvalidTypeForm=false
#!/usr/bin/env python3
"""
PAIR-LRO Verification & Certification Replay Script
===================================================
Verifies all rigorous constants, certified interval bounds (mpmath.iv, 50 dps),
exact algebraic identities, small-cluster reflection positivity / Gaussian domination
checks, and the standard real-hopping RP obstruction.

Claim labels follow the program standard:
  [CERTIFIED] : Rigorous 50-digit interval arithmetic via mpmath.iv
  [EXACT]     : Exact algebraic matrix check on finite Fock spaces
  [NUMERICAL] : Exact diagonalization (ED) on small finite clusters
"""

from __future__ import annotations

import math
import time
import numpy as np
import mpmath as mp
from mpmath import iv

mp.mp.dps = 50
iv.dps = 50


def iv_lo(x) -> mp.mpf:
    return mp.mpf(x.a)


def iv_hi(x) -> mp.mpf:
    return mp.mpf(x.b)


def fmt_floor_sci(x: mp.mpf, sig: int = 5) -> tuple[str, mp.mpf]:
    """Strict directed rounding DOWN (towards 0 / -infinity) of positive x to `sig` significant digits."""
    assert x > 0
    exp = int(mp.floor(mp.log10(x)))
    scaled = x * (mp.mpf(10) ** (sig - 1 - exp))
    mant_int = int(mp.floor(scaled))
    if mant_int >= 10**sig:
        mant_int //= 10
        exp += 1
    s = str(mant_int)
    val_floor = mp.mpf(mant_int) * (mp.mpf(10) ** (exp - sig + 1))
    assert 0 < val_floor <= x
    return f"{s[0]}.{s[1:]}e{exp}", val_floor


def fmt_iv_directed(x_iv, dec: int = 9) -> str:
    """Strict outward directed rounding of an interval [a, b] to `dec` decimal places."""
    scale = mp.mpf(10) ** dec
    lo_int = int(mp.floor(iv_lo(x_iv) * scale))
    hi_int = int(mp.ceil(iv_hi(x_iv) * scale))
    lo_s = f"{lo_int // (10**dec)}.{lo_int % (10**dec):0{dec}d}"
    hi_s = f"{hi_int // (10**dec)}.{hi_int % (10**dec):0{dec}d}"
    return f"[{lo_s}, {hi_s}]"

def section(title: str) -> None:
    print(f"\n{'=' * 72}\n{title}\n{'=' * 72}")


# ---------------------------------------------------------------------------
# 1. [CERTIFIED] Exact Watson integral & 3D box enclosures for W_3, J_sq, I_sin
# ---------------------------------------------------------------------------
section("1. [CERTIFIED] 3D Brillouin-zone integrals (I_W, I_3, W_3(alpha), J_sq(alpha), I_sin)")
t_start = time.time()

I_W_iv = (
    (iv.sqrt(6) / (32 * iv.pi**3))
    * iv.gamma(iv.mpf(1) / 24)
    * iv.gamma(iv.mpf(5) / 24)
    * iv.gamma(iv.mpf(7) / 24)
    * iv.gamma(iv.mpf(11) / 24)
)
print(f"I_W (Watson exact)       in {fmt_iv_directed(I_W_iv, 14)}")


def certify_JW_box(alpha_num: int, alpha_den: int, M: int = 64):
    """
    Rigorous interval enclosure of
      W_3(alpha)   = (1/(2pi)^3) int (alpha - u)_+ / (3*(1+u)) dp
      J_sq(alpha)  = (1/(2pi)^3) int (alpha - u)_+^2 / (1+u) dp
    where u = (cos p1 + cos p2 + cos p3)/3.
    Uses singularity subtraction via exact I_W:
      W_3(alpha)  = ((alpha+1)*I_W - int (u - alpha)_+ / (1+u)) / 3
      J_sq(alpha) = (alpha+1)^2 * I_W - 2*alpha - 1 - int g_2(u; alpha)
    where g_2(u; alpha) = ((u - alpha)_+)^2 / (1+u) is bounded and monotone increasing on (-1, 1].
    Symmetry 0 <= p1 <= p2 <= p3 <= pi reduces the M^3 boxes on [0, pi]^3 by 6x.
    """
    alpha = iv.mpf(alpha_num) / alpha_den
    h = iv.pi / M
    c_pts = [iv.cos(iv.mpf(i) * h) for i in range(M + 1)]
    c_box = [iv.mpf([iv_lo(c_pts[i + 1]), iv_hi(c_pts[i])]) for i in range(M)]

    g1_sum = iv.mpf(0)
    g2_sum = iv.mpf(0)
    for i in range(M):
        ci = c_box[i]
        for j in range(i, M):
            cj = c_box[j]
            for k in range(j, M):
                w = 1 if (i == j == k) else (3 if (i == j or j == k) else 6)
                u_iv = (ci + cj + c_box[k]) / 3
                u_lo, u_hi = iv_lo(u_iv), iv_hi(u_iv)
                # g1(u) = max(u - alpha, 0) / (1 + u), monotone increasing on (-1, 1]
                # g2(u) = max(u - alpha, 0)^2 / (1 + u), monotone increasing on (-1, 1]
                a_val = iv_lo(alpha)
                if u_hi <= a_val:
                    g1_iv = iv.mpf(0)
                    g2_iv = iv.mpf(0)
                else:
                    if u_lo <= a_val:
                        g1_lo = mp.mpf(0)
                        g2_lo = mp.mpf(0)
                    else:
                        d_lo = iv.mpf(u_lo) - alpha
                        one_p_lo = 1 + iv.mpf(u_lo)
                        g1_lo = iv_lo(d_lo / one_p_lo)
                        g2_lo = iv_lo(d_lo * d_lo / one_p_lo)
                    d_hi = iv.mpf(u_hi) - alpha
                    one_p_hi = 1 + iv.mpf(u_hi)
                    g1_hi = iv_hi(d_hi / one_p_hi)
                    g2_hi = iv_hi(d_hi * d_hi / one_p_hi)
                    g1_iv = iv.mpf([g1_lo, g1_hi])
                    g2_iv = iv.mpf([g2_lo, g2_hi])
                g1_sum += w * g1_iv
                g2_sum += w * g2_iv

    inv_M3 = iv.mpf(1) / (M**3)
    g1_int = g1_sum * inv_M3
    g2_int = g2_sum * inv_M3
    W3_iv = ((alpha + 1) * I_W_iv - g1_int) / 3
    Jsq_iv = (alpha + 1) ** 2 * I_W_iv - 2 * alpha - 1 - g2_int
    return W3_iv, Jsq_iv


def certify_Isin_box(M: int = 48):
    """
    Rigorous interval enclosure of I_sin = (1/(2pi)^3) int sqrt(sin^2 p1 + sin^2 p2 + sin^2 p3) dp.
    On [0, pi]^3, sin(p) increases on [0, pi/2] and decreases on [pi/2, pi], and is symmetric
    around pi/2, so [0, pi/2]^3 suffices (with 0 <= p1 <= p2 <= p3 <= pi/2 symmetry).
    """
    h = iv.pi / (2 * M)
    s_pts = [iv.sin(iv.mpf(i) * h) for i in range(M + 1)]
    s2_box = [iv.mpf([iv_lo(s_pts[i]) ** 2, iv_hi(s_pts[i + 1]) ** 2]) for i in range(M)]
    acc = iv.mpf(0)
    for i in range(M):
        si = s2_box[i]
        for j in range(i, M):
            sj = s2_box[j]
            for k in range(j, M):
                w = 1 if (i == j == k) else (3 if (i == j or j == k) else 6)
                acc += w * iv.sqrt(si + sj + s2_box[k])
    return acc / (M**3)


W3_0, Jsq_0 = certify_JW_box(0, 1, M=64)
I3_iv = iv.sqrt(Jsq_0)
W3_half, Jsq_half = certify_JW_box(1, 2, M=64)
W3_3q, Jsq_3q = certify_JW_box(3, 4, M=64)
W3_1 = (2 * I_W_iv - 1) / 3
Jsq_1 = 4 * I_W_iv - 3
I_sin_iv = certify_Isin_box(M=48)

print(f"I_3 = sqrt(J_sq(0))      in {fmt_iv_directed(I3_iv, 9)}  (Koma 2022: 0.68...)")
print(f"W_3(0)                   in {fmt_iv_directed(W3_0, 9)}")
print(f"W_3(1/2)                 in {fmt_iv_directed(W3_half, 9)}")
print(f"W_3(3/4)                 in {fmt_iv_directed(W3_3q, 9)}")
print(f"W_3(1) = (2*I_W - 1)/3   in {fmt_iv_directed(W3_1, 11)}")
print(f"J_sq(0)                  in {fmt_iv_directed(Jsq_0, 9)}")
print(f"J_sq(1/2)                in {fmt_iv_directed(Jsq_half, 9)}")
print(f"J_sq(3/4)                in {fmt_iv_directed(Jsq_3q, 9)}")
print(f"J_sq(1) = 4*I_W - 3      in {fmt_iv_directed(Jsq_1, 11)}")
print(f"I_sin                    in {fmt_iv_directed(I_sin_iv, 9)}  (2*I_sin < sqrt(6)={mp.nstr(mp.sqrt(6), 8)})")
assert iv_hi(2 * I_sin_iv) < mp.sqrt(6)


# ---------------------------------------------------------------------------
# 2. [CERTIFIED] Exact 2-site dimer free energy & rigorous T_0/g certification
# ---------------------------------------------------------------------------
section("2. [CERTIFIED] Full-Fock-space T_0/g bounds (t=0 and pi-flux kappa != 0)")


def H2_iv(r_iv):
    """Certified interval enclosure of binary entropy H2(r) = -r ln r - (1-r) ln(1-r) on [r0, r1]."""
    r0, r1 = iv_lo(r_iv), iv_hi(r_iv)

    def h2_pt(x):
        if x <= 0 or x >= 1:
            return iv.mpf(0)
        xi = iv.mpf(x)
        return -xi * iv.log(xi) - (1 - xi) * iv.log(1 - xi)

    lo = min(iv_lo(h2_pt(r0)), iv_lo(h2_pt(r1)))
    half = mp.mpf("0.5")
    if r0 <= half <= r1:
        hi = iv_hi(iv.log(2))
    else:
        hi = max(iv_hi(h2_pt(r0)), iv_hi(h2_pt(r1)))
    return iv.mpf([lo, hi])


def F_dimer_pair_upper(beta, h_num: int, h_den: int):
    """
    Closed-form certified upper bound on the exact thermodynamic free energy per site
    F_beta / g using the 4-state pair-subspace 2-site AF dimer product ansatz
    supported on P_x x P_y (dim = 4, x in Lambda_even, y in Lambda_odd) with rational
    staggered trial field h = h_num / h_den:
      h_dimer(h)/g = 0.5*(Gamma_x^{(1)} Gamma_y^{(1)} + Gamma_x^{(2)} Gamma_y^{(2)}) - h*(Gamma_x^{(1)} - Gamma_y^{(1)}).

    Exact 4-state spectrum on P_x x P_y (in units of g):
      E_1 = -1/2 - s,  E_2 = 0,  E_3 = -1/2 + s,  E_4 = +1   where s = sqrt(1/4 + 4*h^2).
    Shift by e_max = 1/2 + s:
      Z_4 = 1 + exp(-beta*e_max) + exp(-2*beta*s) + exp(-beta*(e_max + 1))
      m_1 = <Gamma_x^{(1)}>_dimer = -<Gamma_y^{(1)}>_dimer = (2*h/s) * (1 - exp(-2*beta*s)) / Z_4.

    Crucial properties:
      (1) <H_U>_{rho_trial} = 0 identically for ALL U in R (since P_x V_x = 0).
      (2) <H_hop>_{rho_trial} = 0 identically for ALL kappa in R (since single-electron
          hopping flips site fermion parity between even P_x and odd V_x).
      (3) Because rank(rho_dimer) = 4 on 2 sites, the trial entropy per site satisfies
          s_trial = S(rho_dimer)/2 <= (ln 4)/2 = ln 2 <= ln 2 + H_2(r) = S_upper(r)
          for ALL r in [0, 1]. Consequently, for ANY beta >= beta_0:
            E_beta <= e_trial + (S_upper(r) - s_trial) / beta
                   <= e_trial + (S_upper(r) - s_trial) / beta_0
                    = F_dimer_pair(beta_0, h) + (ln 2 + H_2(r)) / beta_0,
          which certifies ALL 0 < T <= T_0 (all beta >= beta_0) in a single evaluation!
    """
    h = iv.mpf(h_num) / h_den
    s = iv.sqrt(iv.mpf(1) / 4 + 4 * h * h)
    e_max = iv.mpf(1) / 2 + s
    Z4 = iv.mpf(1) + iv.exp(-beta * e_max) + iv.exp(-2 * beta * s) + iv.exp(-beta * (e_max + 1))
    lnZ_site = (beta * e_max + iv.log(Z4)) / 2
    if h_num == 0:
        return -lnZ_site / beta
    m1 = (2 * h / s) * (iv.mpf(1) - iv.exp(-2 * beta * s)) / Z4
    return -lnZ_site / beta + h * m1 - iv.mpf(5) / 4 * m1 * m1


def find_best_h_float(beta_f: float) -> tuple[int, int]:
    """Heuristic float search to pick an integer-rational trial field h = p/1000 for rigorous evaluation."""
    best_F = 1e9
    best_p = 0
    for p in range(0, 1501, 10):
        h = p / 1000.0
        s = math.sqrt(0.25 + 4.0 * h * h)
        e_max = 0.5 + s
        Z4 = 1.0 + math.exp(-beta_f * e_max) + math.exp(-2.0 * beta_f * s) + math.exp(-beta_f * (e_max + 1.0))
        lnZ_site = 0.5 * (beta_f * e_max + math.log(Z4))
        m1 = (2.0 * h / s) * (1.0 - math.exp(-2.0 * beta_f * s)) / Z4 if p > 0 else 0.0
        F = -lnZ_site / beta_f + h * m1 - 1.25 * m1 * m1
        if F < best_F:
            best_F = F
            best_p = p
    return best_p, 1000


def certify_all_T_le_T0(
    T0_num: int,
    T0_den: int,
    U_num: int,
    U_den: int,
    kap_num: int,
    kap_den: int,
    alpha_num: int,
    alpha_den: int,
    Jsq_iv,
    nr: int = 800,
) -> mp.mpf:
    """
    Certify positive ODLRO order-parameter lower bound m_0^2 > 0 for ALL 0 < T <= T_0
    (i.e. all beta in [beta_0, infinity)) uniformly in the volume on the full fermion Fock space.
    Uses:
      - Exact 4-state pair-subspace dimer variational bound: since s_trial <= ln 2 <= S_upper(r),
        E_beta <= F_dimer_pair(beta_0, h) + (ln 2 + H_2(r)) / beta_0 for ALL beta >= beta_0.
      - Exact Gaussian domination b_p <= 1 / (beta * g * E_{p+Q})
      - Exact Falk-Bruch + Cauchy-Schwarz bound:
          sum_{p != Q} (alpha-u)_+ g_p / |Lambda| <= sqrt(J_sq(alpha) * (E_1 + |e_hop|/(6g) + I_W/(9*beta_0^2*g^2)))
      - Exact |e_hop| <= min(2*I_sin*|kap|, 2*sqrt(6)*|kap|*sqrt(r))
      - Exact vacancy partition bound at kap=0, U>=0: r(beta) <= 1/(1 + exp(beta*U/2)) <= 1/(1 + exp(beta_0*U/2))
      - Exact 7-site star bound: E_1 <= (1 - r)/sqrt(3)
    Returns certified lower bound on inf_{0 < T <= T_0} m_0^2(T) (strictly > 0 => certified ODLRO for all 0 < T <= T_0).
    """
    T0 = iv.mpf(T0_num) / T0_den
    beta0 = 1 / T0
    U = iv.mpf(U_num) / U_den
    kap = abs(iv.mpf(kap_num) / kap_den)
    alpha = iv.mpf(alpha_num) / alpha_den

    hp, hd = find_best_h_float(float(T0_den / T0_num))
    F_b_0 = F_dimer_pair_upper(beta0, 0, 1)
    F_b_h = F_dimer_pair_upper(beta0, hp, hd)
    F_up = F_b_h if iv_hi(F_b_h) < iv_hi(F_b_0) else F_b_0

    if kap_num == 0 and U_num >= 0:
        r_max = min(mp.mpf(1), iv_hi(iv.mpf(1) / (iv.mpf(1) + iv.exp(beta0 * U / 2))))
    else:
        r_max = mp.mpf(1)

    e_hop_glob = 2 * I_sin_iv * kap
    worst_m2_lo = mp.mpf("1e9")

    for k in range(nr):
        r0 = r_max * k / nr
        r1 = r_max * (k + 1) / nr
        r_iv = iv.mpf([r0, r1])
        S_up = iv.log(2) + H2_iv(r_iv)

        if kap_num != 0:
            e_hop_r = 2 * iv.sqrt(6) * kap * iv.sqrt(r_iv)
            e_hop_iv = iv.mpf([min(iv_lo(e_hop_glob), iv_lo(e_hop_r)), min(iv_hi(e_hop_glob), iv_hi(e_hop_r))])
        else:
            e_hop_iv = iv.mpf(0)

        c_tot = e_hop_iv / 6 + I_W_iv / (9 * beta0 * beta0)
        E_lo_raw = (-F_up - S_up / beta0 + (U * r_iv) / 2 - e_hop_iv) / 3
        E_hi_iv = (1 - r_iv) / iv.sqrt(3)

        if iv_lo(E_lo_raw) > iv_hi(E_hi_iv):
            continue

        E_lo_a = max(mp.mpf(0), iv_lo(E_lo_raw))
        E_hi_b = max(E_lo_a, iv_hi(E_hi_iv))
        E_star_iv = Jsq_iv / 4 - c_tot
        if E_hi_b <= iv_lo(E_star_iv):
            E_iv = iv.mpf(E_hi_b)
            f_iv = E_iv - iv.sqrt(Jsq_iv * (E_iv + c_tot))
        elif E_lo_a >= iv_hi(E_star_iv):
            E_iv = iv.mpf(E_lo_a)
            f_iv = E_iv - iv.sqrt(Jsq_iv * (E_iv + c_tot))
        else:
            f_iv = -Jsq_iv / 4 - c_tot

        m2_iv = (alpha * (1 - r_iv) + f_iv) / (alpha + 1)
        if iv_lo(m2_iv) < worst_m2_lo:
            worst_m2_lo = iv_lo(m2_iv)

    return worst_m2_lo


# Certified (T0/g, alpha) targets verified by 50-digit interval arithmetic over ALL 0 < T <= T0
certified_targets = [
    ("t=0, U/g=-0.5", -1, 2, 0, 1, 155, 1000, 1, 2, Jsq_half),
    ("t=0, U/g=0.0", 0, 1, 0, 1, 205, 1000, 1, 2, Jsq_half),
    ("t=0, U/g=1.0", 1, 1, 0, 1, 435, 1000, 3, 4, Jsq_3q),
    ("t=0, U/g=4.0", 4, 1, 0, 1, 900, 1000, 3, 4, Jsq_3q),
    ("t=0, U/g=10.0", 10, 1, 0, 1, 1140, 1000, 1, 1, Jsq_1),
    ("t=0, U/g=100.0", 100, 1, 0, 1, 1170, 1000, 1, 1, Jsq_1),
    ("pi-flux U=0, kap=0.02", 0, 1, 2, 100, 160, 1000, 1, 2, Jsq_half),
    ("pi-flux U=0, kap=0.05", 0, 1, 5, 100, 85, 1000, 1, 2, Jsq_half),
    ("pi-flux U=4, kap=0.05", 4, 1, 5, 100, 400, 1000, 1, 2, Jsq_half),
    ("pi-flux U=4, kap=0.10", 4, 1, 10, 100, 290, 1000, 1, 2, Jsq_half),
    ("pi-flux U=4, kap=0.20", 4, 1, 20, 100, 45, 1000, 1, 2, Jsq_half),
    ("pi-flux U=10, kap=0.10", 10, 1, 10, 100, 590, 1000, 1, 2, Jsq_half),
    ("pi-flux U=10, kap=0.20", 10, 1, 20, 100, 410, 1000, 1, 2, Jsq_half),
    ("pi-flux U=10, kap=0.35", 10, 1, 35, 100, 80, 1000, 1, 2, Jsq_half),
    ("pi-flux U=20, kap=0.50", 20, 1, 50, 100, 210, 1000, 1, 2, Jsq_half),
]

certified_results: dict[str, tuple[str, mp.mpf, mp.mpf]] = {}
for label, Un, Ud, kn, kd, Tn, Td, an, ad, Jsq_b in certified_targets:
    m2_lo = certify_all_T_le_T0(Tn, Td, Un, Ud, kn, kd, an, ad, Jsq_b)
    floor_str, floor_val = fmt_floor_sci(m2_lo, 5)
    certified_results[label] = (floor_str, floor_val, m2_lo)
    status = "CERTIFIED" if m2_lo > 0 else "FAIL"
    print(f"  [{status}] {label:24s} | alpha={an}/{ad} | 0 < T/g <= {Tn/Td:.4f} | inf m_0^2 >= {floor_str}")
    assert m2_lo > 0 and floor_val <= m2_lo, f"Certification failed for {label}: m2_lo={m2_lo}"


# ---------------------------------------------------------------------------
# 3. [EXACT] Algebraic operator identities on finite fermion Fock spaces
# ---------------------------------------------------------------------------
section("3. [EXACT] Algebraic identities (projector, [H_U, Gamma]=0, double commutators, dimer ED, star bound)")


def make_fermions(n_sites: int):
    n_modes = 2 * n_sites
    dim = 1 << n_modes
    c = []
    for m in range(n_modes):
        mat = np.zeros((dim, dim), dtype=complex)
        for st in range(dim):
            if (st >> m) & 1:
                sgn = -1.0 if ((st & ((1 << m) - 1)).bit_count() % 2 == 1) else 1.0
                mat[st ^ (1 << m), st] = sgn
        c.append(mat)
    cd = [m.conj().T for m in c]

    def a(x, s):
        return c[2 * x + s]

    def ad(x, s):
        return cd[2 * x + s]

    return dim, a, ad


def comm(A, B):
    return A @ B - B @ A


dim2, a2, ad2 = make_fermions(2)
I2 = np.eye(dim2, dtype=complex)
G1 = [ad2(x, 0) @ ad2(x, 1) + a2(x, 1) @ a2(x, 0) for x in range(2)]
G2 = [1j * (ad2(x, 0) @ ad2(x, 1) - a2(x, 1) @ a2(x, 0)) for x in range(2)]
G3 = [I2 - ad2(x, 0) @ a2(x, 0) - ad2(x, 1) @ a2(x, 1) for x in range(2)]
n_up = [ad2(x, 0) @ a2(x, 0) for x in range(2)]
n_dn = [ad2(x, 1) @ a2(x, 1) for x in range(2)]
PV = [n_up[x] + n_dn[x] - 2.0 * (n_up[x] @ n_dn[x]) for x in range(2)]
PP = [I2 - PV[x] for x in range(2)]

# (a) [Gamma_x^(1)]^2 = [Gamma_x^(2)]^2 = P_{P_x}, H_U identity, and Koma square-completion identity
for x in range(2):
    err1 = np.max(np.abs(G1[x] @ G1[x] - PP[x]))
    err2 = np.max(np.abs(G2[x] @ G2[x] - PP[x]))
    HU_op = -(n_up[x] @ n_dn[x] - 0.5 * (n_up[x] + n_dn[x]))
    errU = np.max(np.abs(HU_op - 0.5 * (I2 - G1[x] @ G1[x])))
    assert err1 < 1e-12 and err2 < 1e-12 and errU < 1e-12

sq_comp = 0.25 * ((G1[0] + G1[1]) @ (G1[0] + G1[1]) - (G2[0] - G2[1]) @ (G2[0] - G2[1]))
pair_af_bond = 0.5 * (G1[0] @ G1[1] + G2[0] @ G2[1])
assert np.max(np.abs(sq_comp - pair_af_bond)) < 1e-12

# (b) [H_U, Gamma_y^(a)] = 0 for a in {1, 2, 3}
for x in range(2):
    for y in range(2):
        assert np.max(np.abs(comm(PV[x], G1[y]))) < 1e-12
        assert np.max(np.abs(comm(PV[x], G2[y]))) < 1e-12
        assert np.max(np.abs(comm(PV[x], G3[y]))) < 1e-12

# (c) Complete on-site AND cross-site double commutators on 16-state 2-site Fock space (Lemma 3.1):
# Pair bond P_{0,1} = 0.5*(G1_0 G1_1 + G2_0 G2_1)
dc_pair_onsite = comm(G1[0], comm(pair_af_bond, G1[0])) + comm(G1[1], comm(pair_af_bond, G1[1]))
dc_pair_cross_11 = comm(G1[0], comm(G1[0] @ G1[1], G1[1]))
dc_pair_cross_22 = comm(G1[0], comm(G2[0] @ G2[1], G1[1]))
dc_pair_cross = comm(G1[0], comm(pair_af_bond, G1[1])) + comm(G1[1], comm(pair_af_bond, G1[0]))
assert np.max(np.abs(dc_pair_onsite - (-4.0 * G2[0] @ G2[1]))) < 1e-12
assert np.max(np.abs(dc_pair_cross_11)) < 1e-12
assert np.max(np.abs(dc_pair_cross_22 - (4.0 * G3[0] @ G3[1]))) < 1e-12
assert np.max(np.abs(dc_pair_cross - (4.0 * G3[0] @ G3[1]))) < 1e-12

# Hopping bond K_{0,1} = sum_s i*(ad_0 a_1 - ad_1 a_0) (unit kappa = 1)
K_hop = sum(1j * (ad2(0, s) @ a2(1, s) - ad2(1, s) @ a2(0, s)) for s in (0, 1))
assert np.max(np.abs(PP[0] @ K_hop @ PP[0])) < 1e-12 and np.max(np.abs(PV[0] @ K_hop @ PV[0])) < 1e-12
assert np.max(np.abs(G1[0] @ K_hop @ G1[0])) < 1e-12 and np.max(np.abs(G1[1] @ K_hop @ G1[1])) < 1e-12
dc_hop_onsite = comm(G1[0], comm(K_hop, G1[0])) + comm(G1[1], comm(K_hop, G1[1]))
dc_hop_cross = comm(G1[0], comm(K_hop, G1[1])) + comm(G1[1], comm(K_hop, G1[0]))
dc_hop_cross_formula = 2.0 * (
    G1[0] @ K_hop @ G1[1] + G1[1] @ K_hop @ G1[0] - G1[0] @ G1[1] @ K_hop - K_hop @ G1[0] @ G1[1]
)
assert np.max(np.abs(dc_hop_onsite - (-2.0 * K_hop))) < 1e-12
assert np.max(np.abs(dc_hop_cross - dc_hop_cross_formula)) < 1e-12
assert abs(float(np.linalg.norm(dc_hop_cross, "fro")) - 8.0) < 1e-12

# Repulsive density bond h_V = G3_0 G3_1
h_V = G3[0] @ G3[1]
dc_V_onsite = comm(G1[0], comm(h_V, G1[0])) + comm(G1[1], comm(h_V, G1[1]))
dc_V_cross = comm(G1[0], comm(h_V, G1[1])) + comm(G1[1], comm(h_V, G1[0]))
assert np.max(np.abs(dc_V_onsite - (-8.0 * h_V))) < 1e-12
assert np.max(np.abs(dc_V_cross - (8.0 * G2[0] @ G2[1]))) < 1e-12

# (d) Exact 2-site pair-subspace ED check of F_dimer_pair_upper (Lemma 4.2, E_2 = 0):
pair_idx_2 = [st for st in range(dim2) if np.real(PP[0][st, st] * PP[1][st, st]) > 0.5]
assert len(pair_idx_2) == 4
P_sub_op = pair_af_bond[np.ix_(pair_idx_2, pair_idx_2)]
G1_0_sub = G1[0][np.ix_(pair_idx_2, pair_idx_2)]
G1_1_sub = G1[1][np.ix_(pair_idx_2, pair_idx_2)]
for h_num_t, h_den_t in [(0, 1), (1, 2), (3, 4), (3, 2)]:
    h_t = h_num_t / h_den_t
    s_t = math.sqrt(0.25 + 4.0 * h_t * h_t)
    H_dim_sub = P_sub_op - h_t * (G1_0_sub - G1_1_sub)
    evals_dim, evecs_dim = np.linalg.eigh(H_dim_sub)
    expected_evals = np.sort([-0.5 - s_t, 0.0, -0.5 + s_t, 1.0])
    assert np.max(np.abs(evals_dim - expected_evals)) < 1e-12
    for beta_t in [1.0, 2.5, 5.0]:
        w_d = np.exp(-beta_t * evals_dim)
        rho_d = (evecs_dim * (w_d / np.sum(w_d))) @ evecs_dim.conj().T
        mx_d = float(np.real(np.trace(rho_d @ G1_0_sub)))
        my_d = float(np.real(np.trace(rho_d @ G1_1_sub)))
        assert abs(mx_d + my_d) < 1e-12
        e_intra_site = 0.5 * float(np.real(np.trace(rho_d @ P_sub_op)))
        e_inter_site = 2.5 * 0.5 * mx_d * my_d  # 5/2 inter-dimer bonds/site * (g/2) * m_x * m_y = -1.25 * m_1^2
        s_site = -0.5 * float(np.sum((w_d / np.sum(w_d)) * np.log(w_d / np.sum(w_d))))
        F_ed_site = e_intra_site + e_inter_site - s_site / beta_t
        F_iv_site = F_dimer_pair_upper(iv.mpf(str(beta_t)), h_num_t, h_den_t)
        assert abs(F_ed_site - float(iv_lo(F_iv_site))) < 1e-13

# (e) Exact 7-site star bound || 0.5 * sum_{j=1}^6 (G1_0 G1_j + G2_0 G2_j) || = 2*sqrt(3)
# Verified BOTH in the reduced 14x14 (S=3) sector AND on the full 2^7 = 128-state pair subspace P^{\otimes 7}
Sp0 = np.array([[0.0, 1.0], [0.0, 0.0]])
Sm0 = Sp0.T
S3_p = np.zeros((7, 7))
for m_idx in range(6):
    m_val = 3 - m_idx
    S3_p[m_idx, m_idx + 1] = math.sqrt(3 * 4 - m_val * (m_val - 1))
S3_m = S3_p.T
star_op_14 = np.kron(Sp0, S3_m) + np.kron(Sm0, S3_p)
star_norm_14 = float(np.linalg.eigvalsh(star_op_14)[-1])
assert abs(star_norm_14 - 2.0 * math.sqrt(3.0)) < 1e-12

sx_2 = np.array([[0.0, 1.0], [1.0, 0.0]])
sy_2 = np.array([[0.0, -1.0j], [1.0j, 0.0]])
I_2 = np.eye(2)
star_op_128 = np.zeros((128, 128), dtype=complex)
for nb in range(1, 7):
    ops_x = [sx_2 if idx in (0, nb) else I_2 for idx in range(7)]
    ops_y = [sy_2 if idx in (0, nb) else I_2 for idx in range(7)]
    term_x, term_y = ops_x[0], ops_y[0]
    for idx in range(1, 7):
        term_x = np.kron(term_x, ops_x[idx])
        term_y = np.kron(term_y, ops_y[idx])
    star_op_128 += -0.5 * (term_x + term_y)
star_norm_128 = float(np.linalg.eigvalsh(star_op_128)[-1])
assert abs(star_norm_128 - 2.0 * math.sqrt(3.0)) < 1e-12


# (f) Exact 2-site verification of Lemma 3.2 Eqs. (3.11)-(3.13) vacancy-hopping decomposition:
B2_hop = [[1j * a2(1, s) for s in (0, 1)], [-1j * a2(0, s) for s in (0, 1)]]
for x in range(2):
    for s in (0, 1):
        assert float(np.linalg.norm(ad2(x, s) - (PV[x] @ ad2(x, s) + ad2(x, s) @ PV[x]), "fro")) < 1e-12
# Eq. (3.11): H_hop = T_1 + T_2
T1_hop = sum(PV[x] @ ad2(x, s) @ B2_hop[x][s] for x in range(2) for s in (0, 1))
T2_hop = sum(ad2(x, s) @ PV[x] @ B2_hop[x][s] for x in range(2) for s in (0, 1))
assert float(np.linalg.norm(K_hop - (T1_hop + T2_hop), "fro")) < 1e-12
# Eq. (3.12): T_1 and T_2^\dagger forms
T1_rhs = sum((a2(x, s) @ PV[x]).conj().T @ B2_hop[x][s] for x in range(2) for s in (0, 1))
T2_dag_1 = sum(B2_hop[x][s].conj().T @ PV[x] @ a2(x, s) for x in range(2) for s in (0, 1))
T2_dag_2 = -sum((PV[x] @ a2(x, s)) @ B2_hop[x][s].conj().T for x in range(2) for s in (0, 1))
T2_dag_3 = -sum((ad2(x, s) @ PV[x]).conj().T @ B2_hop[x][s].conj().T for x in range(2) for s in (0, 1))
assert float(np.linalg.norm(T1_hop - T1_rhs, "fro")) < 1e-12
assert float(np.linalg.norm(T2_hop.conj().T - T2_dag_1, "fro")) < 1e-12
assert float(np.linalg.norm(T2_hop.conj().T - T2_dag_2, "fro")) < 1e-12
assert float(np.linalg.norm(T2_hop.conj().T - T2_dag_3, "fro")) < 1e-12
# Eq. (3.13) operator identities and Cauchy-Schwarz trace inequality on 2-site Fock space
for x in range(2):
    sum_p_ad_a_p = sum(PV[x] @ ad2(x, s) @ a2(x, s) @ PV[x] for s in (0, 1))
    sum_p_a_ad_p = sum(PV[x] @ a2(x, s) @ ad2(x, s) @ PV[x] for s in (0, 1))
    assert float(np.linalg.norm(sum_p_ad_a_p - PV[x], "fro")) < 1e-12
    assert float(np.linalg.norm(sum_p_a_ad_p - PV[x], "fro")) < 1e-12
H2_test = pair_af_bond + 0.7 * K_hop + 0.3 * (G1[0] - G1[1])
ev2_t, evec2_t = np.linalg.eigh(H2_test)
w2_t = np.exp(-1.5 * ev2_t)
rho2_t = (evec2_t * (w2_t / np.sum(w2_t))) @ evec2_t.conj().T
lhs_313 = abs(complex(np.trace(rho2_t @ K_hop)))
f1_a = math.sqrt(sum(float(np.real(np.trace(rho2_t @ PV[x] @ ad2(x, s) @ a2(x, s) @ PV[x]))) for x in range(2) for s in (0, 1)))
f2_a = math.sqrt(sum(float(np.real(np.trace(rho2_t @ B2_hop[x][s].conj().T @ B2_hop[x][s]))) for x in range(2) for s in (0, 1)))
f1_b = math.sqrt(sum(float(np.real(np.trace(rho2_t @ PV[x] @ a2(x, s) @ ad2(x, s) @ PV[x]))) for x in range(2) for s in (0, 1)))
f2_b = math.sqrt(sum(float(np.real(np.trace(rho2_t @ B2_hop[x][s] @ B2_hop[x][s].conj().T))) for x in range(2) for s in (0, 1)))
assert lhs_313 <= f1_a * f2_a + f1_b * f2_b + 1e-12
print(
    "  [EXACT] All algebraic identities (3.1)-(3.13), on/cross-site double commutators, 4-state dimer ED, and 128-state star norm (2*sqrt(3)) verified."
)

# ---------------------------------------------------------------------------
# 4. [NUMERICAL] 1D 4-site PBC: Gaussian domination & factor-of-2 correction check
# ---------------------------------------------------------------------------
section("4. [NUMERICAL] 4-site ring (dim=256): Koma U=0 AF & U=1 FM Duhamel factor-of-2 check + p-sum identity")

dim4, a4, ad4 = make_fermions(4)
I4 = np.eye(dim4, dtype=complex)
G1_4 = [ad4(x, 0) @ ad4(x, 1) + a4(x, 1) @ a4(x, 0) for x in range(4)]
G2_4 = [1j * (ad4(x, 0) @ ad4(x, 1) - a4(x, 1) @ a4(x, 0)) for x in range(4)]
PV_4 = [
    (ad4(x, 0) @ a4(x, 0)) + (ad4(x, 1) @ a4(x, 1)) - 2.0 * (ad4(x, 0) @ a4(x, 0) @ ad4(x, 1) @ a4(x, 1))
    for x in range(4)
]

# Verify momentum-space Cauchy-Schwarz vacancy-hopping bound (Lemma 3.2) on 4-site ring (d=1):
# |<H_hop>| / N <= 2 * sqrt(2*d) * |kap| * sqrt(r) with d=1 => 2 * sqrt(2) * |kap| * sqrt(r)
H4_hop_koma = sum(
    1j * 0.1 * (1.0 if x < 3 else -1.0) * (ad4(x, s) @ a4((x + 1) % 4, s) - ad4((x + 1) % 4, s) @ a4(x, s))
    for x in range(4)
    for s in (0, 1)
)
A_cs = np.zeros((dim4, dim4), dtype=complex)
B_cs = np.zeros((dim4, dim4), dtype=complex)
for x in range(4):
    s_x = 1.0 if x < 3 else -1.0
    s_xm = 1.0 if ((x - 1) % 4) < 3 else -1.0
    for s in (0, 1):
        Bx = 1j * 0.1 * (s_x * a4((x + 1) % 4, s) - s_xm * a4((x - 1) % 4, s))
        A_cs += Bx.conj().T @ Bx
        B_cs += Bx @ Bx.conj().T
assert np.max(np.abs((A_cs + B_cs) - (4.0 * 1.0 * (0.1**2) * 4.0 * I4))) < 1e-12

# Verify exact p-averaged double-commutator identity on 4-site ring:
# (1/4) sum_p [Gamma_{-p}^(1), [H, Gamma_p^(1)]] == (1/4) sum_x [Gamma_x^(1), [H, Gamma_x^(1)]]
H4_pair_AF = sum(0.5 * (G1_4[x] @ G1_4[(x + 1) % 4] + G2_4[x] @ G2_4[(x + 1) % 4]) for x in range(4))
H4_pair_FM = -H4_pair_AF
H4_U = sum(0.5 * (I4 - G1_4[x] @ G1_4[x]) for x in range(4))
H4_hop_FM = sum(
    -0.1 * (1.0 if x < 3 else -1.0) * (ad4(x, s) @ a4((x + 1) % 4, s) + ad4((x + 1) % 4, s) @ a4(x, s))
    for x in range(4)
    for s in (0, 1)
)

Gp_list = [sum(np.exp(-1j * (2.0 * math.pi * k_idx / 4.0) * x) * G1_4[x] for x in range(4)) / 2.0 for k_idx in range(4)]
for H_test in [H4_pair_AF + H4_hop_koma, H4_pair_FM + 1.0 * H4_U + H4_hop_FM]:
    dc_p_avg = sum(comm(Gp.conj().T, comm(H_test, Gp)) for Gp in Gp_list) / 4.0
    dc_x_avg = sum(comm(G1_4[x], comm(H_test, G1_4[x])) for x in range(4)) / 4.0
    assert np.max(np.abs(dc_p_avg - dc_x_avg)) < 1e-12


def check_4site_duhamel(H_mat: np.ndarray, Q_shift: float, label_str: str):
    evals, evecs = np.linalg.eigh(H_mat)
    evals = evals - evals[0]
    for beta_val in [1.0, 2.0, 5.0]:
        w = np.exp(-beta_val * evals)
        Z = float(np.sum(w))
        dE = evals[None, :] - evals[:, None]
        dw = w[:, None] - w[None, :]
        K = np.where(np.abs(dE) > 1e-10, dw / np.where(np.abs(dE) > 1e-10, beta_val * dE, 1.0), w[:, None])
        max_true_ratio = 0.0
        max_koma_ratio = 0.0
        for k_idx in range(4):
            p = 2.0 * math.pi * k_idx / 4.0
            EpQ = 1.0 - math.cos(p + Q_shift)
            if EpQ < 1e-9:
                continue
            Gp = Gp_list[k_idx]
            Ad = evecs.conj().T @ Gp.conj().T @ evecs
            Be = evecs.conj().T @ Gp @ evecs
            b_p = float(np.sum(K * Ad * Be.T).real) / Z
            c_p = float(np.real(np.sum(w * np.diag(evecs.conj().T @ comm(Gp.conj().T, comm(H_mat, Gp)) @ evecs)))) / Z
            assert c_p >= -1e-12
            true_bound = 1.0 / (beta_val * 1.0 * EpQ)
            koma_bound = 1.0 / (2.0 * beta_val * 1.0 * EpQ)
            max_true_ratio = max(max_true_ratio, b_p / true_bound)
            max_koma_ratio = max(max_koma_ratio, b_p / koma_bound)
        print(
            f"  [{label_str}] beta={beta_val:3.1f}: max b_p/(1/(beta*g*E)) = {max_true_ratio:.6f} <= 1 | vs Koma 1/(2*beta*g*E) = {max_koma_ratio:.6f}"
        )
        assert max_true_ratio <= 1.0 + 1e-9
    return evals, evecs


# (i) Koma's exact U=0 AF model (Q = pi) at kap = 0.0 and kap = 0.1
evals4_koma0, evecs4_koma0 = check_4site_duhamel(H4_pair_AF, math.pi, "Koma U=0 AF, kap=0.0")
evals4_koma, evecs4_koma = check_4site_duhamel(H4_pair_AF + H4_hop_koma, math.pi, "Koma U=0 AF, kap=0.1")
# (ii) Full FM + U=1 + pi-flux model (Q = 0) with kappa = 0.1
evals4, evecs4 = check_4site_duhamel(H4_pair_FM + 1.0 * H4_U + H4_hop_FM, 0.0, "FM + U=1,    kap=0.1")


# ---------------------------------------------------------------------------
# 5. [NUMERICAL] Obstruction witness: Majorana RP fails for 0-flux real hopping
# ---------------------------------------------------------------------------
section("5. [NUMERICAL] 2x2 Plaquette (14x14 Gram matrix): Majorana RP holds for pi-flux, FAILS for 0-flux")

# 4 sites: 0=(0,0) [even], 1=(0,1) [odd], 2=(1,0) [odd], 3=(1,1) [even]
# Cut across x^(1): Lambda_- = {0,1}, Lambda_+ = {2,3}; reflection theta swaps 0<->2 and 1<->3.
n_modes_p = 8
dim_p = 1 << n_modes_p
cp, cdp = [], []
for m in range(n_modes_p):
    mat = np.zeros((dim_p, dim_p), dtype=complex)
    for st in range(dim_p):
        if (st >> m) & 1:
            sgn = -1.0 if ((st & ((1 << m) - 1)).bit_count() % 2 == 1) else 1.0
            mat[st ^ (1 << m), st] = sgn
    cp.append(mat)
    cdp.append(mat.conj().T)
ap = lambda x, s: cp[2 * x + s]
adp = lambda x, s: cdp[2 * x + s]
np_op = lambda x, s: adp(x, s) @ ap(x, s)

Ip = np.eye(dim_p, dtype=complex)
G1p = [adp(x, 0) @ adp(x, 1) + ap(x, 1) @ ap(x, 0) for x in range(4)]
G2p = [1j * (adp(x, 0) @ adp(x, 1) - ap(x, 1) @ ap(x, 0)) for x in range(4)]
xip = [[ap(x, s) + adp(x, s) for s in (0, 1)] for x in range(4)]
etap = [[1j * (adp(x, s) - ap(x, s)) for s in (0, 1)] for x in range(4)]

# Antilinear reflection theta: theta(xi_{x,s}) = xi_{th(x),s}, theta(eta_{x,s}) = -eta_{th(x),s}
th_site = {0: 2, 1: 3, 2: 0, 3: 1}
perm_modes = [2 * th_site[x] + s for x in range(4) for s in range(2)]
S_th = np.zeros((dim_p, dim_p), dtype=complex)
for st in range(dim_p):
    vec = np.zeros(dim_p, dtype=complex)
    vec[0] = 1.0
    for m in reversed([m for m in range(n_modes_p) if (st >> m) & 1]):
        vec = cdp[perm_modes[m]] @ vec
    S_th[:, st] = vec


def theta_map(A: np.ndarray) -> np.ndarray:
    return S_th @ A.conj() @ S_th.T


# Koma's exact bipartite particle-hole unitary U_ph_odd on odd sites {1, 2}
# (flips eta_{x,s} -> -eta_{x,s}, leaves xi_{x,s} unchanged, flips Gamma^(1) -> -Gamma^(1) and Gamma^(3) -> -Gamma^(3))
U_ph_odd = np.eye(dim_p, dtype=complex)
for x in (1, 2):
    for s in (0, 1):
        par = np.eye(dim_p, dtype=complex)
        for y in range(4):
            for sp in (0, 1):
                if y != x or sp != s:
                    par = par @ (Ip - 2.0 * np_op(y, sp))
        U_ph_odd = U_ph_odd @ (par @ (adp(x, s) + ap(x, s)))

# Bipartite pi/2 gauge U_phase_odd on odd sites {1, 2}: relates FM pair + real pi-flux <-> Koma AF pair + imag pi-flux
U_phase_odd = np.eye(dim_p, dtype=complex)
for x in (1, 2):
    for s in (0, 1):
        # exp(i*pi/2 * n) diagonal on occupation states
        U_phase_odd = U_phase_odd @ np.diag(
            [np.exp(1j * 0.5 * math.pi * ((st >> (2 * x + s)) & 1)) for st in range(dim_p)]
        )

ops_minus = [
    Ip,
    G1p[0],
    G2p[0],
    G1p[1],
    G2p[1],
    G1p[0] @ G1p[1],
    G2p[0] @ G2p[1],
    xip[0][0] @ xip[0][1],
    etap[0][0] @ etap[0][1],
    xip[0][0] @ xip[1][0],
    etap[0][0] @ etap[1][0],
    np_op(0, 0),
    np_op(1, 0),
    adp(0, 0) @ ap(1, 0) + adp(1, 0) @ ap(0, 0),
]
assert len(ops_minus) == 14
thetas_minus = [theta_map(op) for op in ops_minus]


def test_rp_plaquette(H_orig: np.ndarray, U_unitary: np.ndarray, beta_val: float) -> tuple[float, float]:
    H_tilde = U_unitary @ H_orig @ U_unitary.conj().T
    evals, evecs = np.linalg.eigh(H_tilde)
    rho = (evecs * np.exp(-beta_val * (evals - evals[0]))) @ evecs.conj().T
    rho /= np.trace(rho).real
    G = np.zeros((len(ops_minus), len(ops_minus)), dtype=complex)
    for i, A in enumerate(ops_minus):
        for j, thB in enumerate(thetas_minus):
            G[i, j] = np.trace(A @ thB @ rho)
    herm_err = float(np.max(np.abs(G - G.conj().T)))
    min_eig = float(np.min(np.linalg.eigvalsh(0.5 * (G + G.conj().T))))
    return herm_err, min_eig


plaq_edges = [(0, 2), (1, 3), (0, 1), (2, 3)]
H_plaq_pair_AF = sum(+0.5 * (G1p[x] @ G1p[y] + G2p[x] @ G2p[y]) for x, y in plaq_edges)
H_plaq_pair_FM = -H_plaq_pair_AF
H_plaq_U = sum(0.5 * (Ip - G1p[x] @ G1p[x]) for x in range(4))

# Koma pi-flux: i*kap*(ad_x a_y - h.c.) with signs (+1, +1, +1, -1) -> flux pi per plaquette
pi_koma_signs = [(0, 2, 1.0), (1, 3, 1.0), (0, 1, 1.0), (2, 3, -1.0)]
H_plaq_pi_koma = sum(
    1j * 0.5 * sgn * (adp(x, s) @ ap(y, s) - adp(y, s) @ ap(x, s)) for x, y, sgn in pi_koma_signs for s in (0, 1)
)
H_plaq_pi_real = U_phase_odd.conj().T @ H_plaq_pi_koma @ U_phase_odd
H_plaq_real0 = sum(-0.5 * (adp(x, s) @ ap(y, s) + adp(y, s) @ ap(x, s)) for x, y in plaq_edges for s in (0, 1))

h_koma, me_koma = test_rp_plaquette(H_plaq_pair_AF + H_plaq_U + H_plaq_pi_koma, U_ph_odd, beta_val=2.0)
h_fm_pi, me_fm_pi = test_rp_plaquette(H_plaq_pair_FM + H_plaq_U + H_plaq_pi_real, U_ph_odd @ U_phase_odd, beta_val=2.0)
h_r0a, me_r0a = test_rp_plaquette(H_plaq_pair_FM + H_plaq_U + H_plaq_real0, U_ph_odd @ U_phase_odd, beta_val=2.0)
h_r0b, me_r0b = test_rp_plaquette(H_plaq_pair_FM + H_plaq_U + H_plaq_real0, U_ph_odd, beta_val=2.0)

print(f"  2x2 AF pair + U + Koma pi-flux (U_ph):           herm_err={h_koma:.2e}, min_eig={me_koma:+.6e} (>= 0: RP HOLDS)")
print(f"  2x2 FM pair + U + real pi-flux (U_ph*U_phase):   herm_err={h_fm_pi:.2e}, min_eig={me_fm_pi:+.6e} (>= 0: RP HOLDS)")
print(f"  2x2 FM pair + U + real 0-flux  (U_ph*U_phase):   herm_err={h_r0a:.2e}, min_eig={me_r0a:+.6e} (< 0: RP FAILS)")
print(f"  2x2 FM pair + U + real 0-flux  (U_ph only):      herm_err={h_r0b:.2e}, min_eig={me_r0b:+.6e} (< 0: RP FAILS)")
assert h_koma < 1e-12 and me_koma >= -1e-12
assert h_fm_pi < 1e-12 and me_fm_pi >= -1e-12
assert me_r0a < -1e-2 and me_r0b < -1e-2

# ---------------------------------------------------------------------------
# 6. [CERTIFIED] Literal self-check of PROOF.md Eqs. (5.5)-(5.9) vs. Tables 1 & 2
# ---------------------------------------------------------------------------
section("6. [CERTIFIED] Literal self-check of PROOF.md Eqs. (5.5)-(5.9) vs. Tables 1 & 2")
# Verify PROOF.md Eq. (3.18) and Section 5.1 definition of J_sq(alpha) = (2pi)^{-3} int (alpha - u_k)_+^2 / (1 + u_k) d^3k
# at alpha = 1 matches J_sq(1) = 4 * I_W - 3 (without any spurious factor of 2):
_, Jsq_1_from_box = certify_JW_box(1, 1, M=16)
assert iv_lo(Jsq_1_from_box) == iv_lo(4 * I_W_iv - 3) and iv_hi(Jsq_1_from_box) == iv_hi(4 * I_W_iv - 3)
grid_1d = (np.arange(32) + 0.5) * (math.pi / 32.0)
c_grid = np.cos(grid_1d)
u_grid = (c_grid[:, None, None] + c_grid[None, :, None] + c_grid[None, None, :]) / 3.0
# Singularity-subtracted definition from PROOF.md Section 5.1:
# J_sq(alpha) = (2pi)^{-3} int [((alpha - u)_+^2 - (1 + alpha)^2) / (1 + u)] d^3k + (1 + alpha)^2 I_W
jsq1_sub_num = float(np.mean(((1.0 - u_grid) ** 2 - 4.0) / (1.0 + u_grid))) + 4.0 * float(iv_lo(I_W_iv))
assert abs(jsq1_sub_num - float(iv_lo(4 * I_W_iv - 3))) < 1e-12

# Exact table rows as displayed in PROOF.md and RESULT.md (Table 1: 6 rows; Table 2: 9 rows)
# Format: (label, U/g_str, |kap|/g_str, T0/g_str, displayed_floor_str, alpha_str)
displayed_table_rows = [
    ("t=0, U/g=-0.5", "-0.5", "0", "0.155", "1.7286e-3", "1/2"),
    ("t=0, U/g=0.0", "0.0", "0", "0.205", "1.6906e-3", "1/2"),
    ("t=0, U/g=1.0", "1.0", "0", "0.435", "2.2130e-3", "3/4"),
    ("t=0, U/g=4.0", "4.0", "0", "0.900", "2.4383e-3", "3/4"),
    ("t=0, U/g=10.0", "10.0", "0", "1.140", "1.1750e-3", "1/1"),
    ("t=0, U/g=100.0", "100.0", "0", "1.170", "1.4858e-3", "1/1"),
    ("pi-flux U=0, kap=0.02", "0.0", "0.02", "0.160", "7.8445e-4", "1/2"),
    ("pi-flux U=0, kap=0.05", "0.0", "0.05", "0.085", "5.8069e-4", "1/2"),
    ("pi-flux U=4, kap=0.05", "4.0", "0.05", "0.400", "1.1635e-3", "1/2"),
    ("pi-flux U=4, kap=0.10", "4.0", "0.10", "0.290", "7.5765e-4", "1/2"),
    ("pi-flux U=4, kap=0.20", "4.0", "0.20", "0.045", "7.8083e-4", "1/2"),
    ("pi-flux U=10, kap=0.10", "10.0", "0.10", "0.590", "7.9024e-4", "1/2"),
    ("pi-flux U=10, kap=0.20", "10.0", "0.20", "0.410", "1.5236e-4", "1/2"),
    ("pi-flux U=10, kap=0.35", "10.0", "0.35", "0.080", "6.5047e-6", "1/2"),
    ("pi-flux U=20, kap=0.50", "20.0", "0.50", "0.210", "4.0079e-4", "1/2"),
]

assert len(displayed_table_rows) == len(certified_targets) == 15

for (label, Un, Ud, kn, kd, Tn, Td, an, ad, Jsq_iv), (
    t_label,
    t_U,
    t_kap,
    t_T0,
    t_floor_str,
    t_alpha,
) in zip(certified_targets, displayed_table_rows):
    assert label == t_label
    assert f"{an}/{ad}" == t_alpha, f"Alpha mismatch on {label}: {an}/{ad} vs {t_alpha}"
    assert abs(Un / Ud - float(t_U)) < 1e-12
    assert abs(kn / kd - float(t_kap)) < 1e-12
    assert abs(Tn / Td - float(t_T0)) < 1e-12

    # Literal evaluation of PROOF.md Eqs. (5.5)-(5.9):
    T0_iv = iv.mpf(Tn) / Td
    beta0_iv = 1 / T0_iv
    U_iv = iv.mpf(Un) / Ud
    kap_iv = iv.mpf(kn) / kd
    alpha_iv = iv.mpf(an) / ad
    hn, hd = find_best_h_float(float(Td / Tn))
    F_b_0 = F_dimer_pair_upper(beta0_iv, 0, 1)
    F_b_h = F_dimer_pair_upper(beta0_iv, hn, hd)
    F_up_iv = F_b_h if iv_hi(F_b_h) < iv_hi(F_b_0) else F_b_0

    if Un >= 0 and kn == 0:
        r_max_iv = iv.mpf(1) / (iv.mpf(1) + iv.exp(beta0_iv * U_iv / 2))
        r_top = min(mp.mpf(1), iv_hi(r_max_iv))
    else:
        r_top = mp.mpf(1)

    e_hop_glob = 2 * I_sin_iv * kap_iv
    nr_check = 800
    literal_min_lo = mp.mpf("1e9")
    for k_sub in range(nr_check):
        r_iv = iv.mpf([r_top * k_sub / nr_check, r_top * (k_sub + 1) / nr_check])
        S_up_iv = iv.log(2) + H2_iv(r_iv)
        if kn != 0:
            e_hop_r = 2 * iv.sqrt(6) * kap_iv * iv.sqrt(r_iv)
            e_hop_up_iv = iv.mpf([min(iv_lo(e_hop_glob), iv_lo(e_hop_r)), min(iv_hi(e_hop_glob), iv_hi(e_hop_r))])
        else:
            e_hop_up_iv = iv.mpf(0)
        # Eq. (5.5) / Eq. (5.8): E_1^lo(r) and E_1^hi(r)
        c_term_iv = e_hop_up_iv / 6 + I_W_iv / (9 * beta0_iv * beta0_iv)
        E1_lo_raw = (-F_up_iv - S_up_iv / beta0_iv + (U_iv * r_iv) / 2 - e_hop_up_iv) / 3
        E1_hi_iv = (1 - r_iv) / iv.sqrt(3)
        if iv_lo(E1_lo_raw) > iv_hi(E1_hi_iv):
            continue
        E_lo_a = max(mp.mpf(0), iv_lo(E1_lo_raw))
        E_hi_b = max(E_lo_a, iv_hi(E1_hi_iv))
        E_star_iv = Jsq_iv / 4 - c_term_iv
        if E_hi_b <= iv_lo(E_star_iv):
            E_eval = iv.mpf(E_hi_b)
            rhs_min = E_eval - iv.sqrt(Jsq_iv * (E_eval + c_term_iv))
        elif E_lo_a >= iv_hi(E_star_iv):
            E_eval = iv.mpf(E_lo_a)
            rhs_min = E_eval - iv.sqrt(Jsq_iv * (E_eval + c_term_iv))
        else:
            rhs_min = -Jsq_iv / 4 - c_term_iv
        m2_sub_iv = (alpha_iv * (1 - r_iv) + rhs_min) / (alpha_iv + 1)
        if iv_lo(m2_sub_iv) < literal_min_lo:
            literal_min_lo = iv_lo(m2_sub_iv)

    floor_str_expected, floor_val_expected, m2_lo_ref = certified_results[label]
    assert literal_min_lo == m2_lo_ref, f"Formula mismatch on {label}: {literal_min_lo} vs {m2_lo_ref}"
    assert t_floor_str == floor_str_expected, f"Table floor string mismatch on {label}: {t_floor_str} vs {floor_str_expected}"
    table_val = mp.mpf(t_floor_str)
    assert 0 < table_val <= literal_min_lo, f"Table value {table_val} exceeds certified lower bound {literal_min_lo}!"

print("  [CERTIFIED] All 15 rows of PROOF.md Eqs. (5.5)-(5.9) & Tables 1-2 match literal evaluation and directed rounding down.")

print(f"\nALL REPLAY CHECKS PASSED in {time.time() - t_start:.2f}s.")
