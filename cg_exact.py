"""
Authoritative Clebsch-Gordan coefficients for the Hg 546.1 nm Zeeman multiplet,
obtained by numerical diagonalisation of the total angular momentum operator
J^2 = J1^2 + J2^2 + 2 J1.J2  in the uncoupled basis |j1 m1, j2 m2>.

This construction is self-validating: the eigenvalues of J^2 must be exactly
J(J+1) with multiplicity 2J+1, which is checked before the coefficients are used.

Emission: upper 6s7s 3S1 (J=1)  ->  lower 6s6p 3P2 (J=2)
The photon carries j2 = 1, q = m_u - m_l.
"""
import numpy as np

np.set_printoptions(suppress=True)


def cg_from_diagonalisation(j1, j2):
    """
    Return dict[(m1, q, M)] = <j1 m1; j2 q | J M> for J = |j1-j2| ... j1+j2,
    built by diagonalising J^2.  Normalisation and phase are fixed by the
    requirement that every fixed-m1 slice is a unitary rotation.
    """
    m1s = [m for m in range(-j1, j1 + 1)]
    q2s = [q for q in range(-j2, j2 + 1)]
    keys = [(m, q) for m in m1s for q in q2s]
    idx = {k: a for a, k in enumerate(keys)}
    n = len(keys)

    # J^2 = J1^2 + J2^2 + 2 (J1z J2z + (J1+J2- + J1-J2+)/2)
    H = np.zeros((n, n))
    for (m, q), a in idx.items():
        H[a, a] = j1 * (j1 + 1) + j2 * (j2 + 1) + 2.0 * m * q
        b = idx.get((m + 1, q - 1))
        if b is not None:
            c = np.sqrt((j1 - m) * (j1 + m + 1)) * np.sqrt((j2 + q) * (j2 - q + 1))
            H[a, b] += c          # (1/2) * 2 * c  from 2 J1.J2
        b = idx.get((m - 1, q + 1))
        if b is not None:
            c = np.sqrt((j1 + m) * (j1 - m + 1)) * np.sqrt((j2 - q) * (j2 + q + 1))
            H[a, b] += c

    assert np.allclose(H, H.T), "Hamiltonian not symmetric"

    ev, evec = np.linalg.eigh(H)

    # ---- validation: spectrum must be J(J+1), multiplicity 2J+1 ----
    Js = list(range(abs(j1 - j2), j1 + j2 + 1))
    expected = []
    for J in Js:
        expected += [J * (J + 1)] * (2 * J + 1)
    expected = np.sort(np.array(expected, dtype=float))
    assert len(expected) == n
    assert np.allclose(ev, expected, atol=1e-9), (
        f"spectrum mismatch\n got {np.round(ev, 5)}\n want {expected}")

    # ---- collect coefficients, grouping eigenvectors by J ----
    # eigh returns an ARBITRARY basis inside each degenerate J block, so we
    # must project onto definite M.  J_z is already diagonal in this basis
    # (J_z = m1 + q), so |J M> = P_J |M> / ||P_J |M>||.
    out = {}
    i = 0
    for J in Js:
        mult = 2 * J + 1
        block = evec[:, i:i + mult]
        i += mult
        for M in range(-J, J + 1):
            # |M> has amplitude 1 on every basis state with m1 + q == M
            amp = np.zeros(n)
            for (m, q), a in idx.items():
                if m + q == M:
                    amp[a] = 1.0
            v = block @ (block.T @ amp)      # projection of |M> onto the J block
            nrm = np.linalg.norm(v)
            if nrm < 1e-9:
                continue
            v = v / nrm
            # fix the overall sign: make the largest component positive
            k = int(np.argmax(np.abs(v)))
            if v[k] < 0:
                v = -v
            for (m, q), a in idx.items():
                if abs(v[a]) > 1e-9:
                    out[(m, q, M)] = v[a]
    return out


def check_orthonormality(out, j1, j2):
    """For each fixed m1, sum_{q,M} |CG|^2 must equal 1."""
    worst = 0.0
    for m in range(-j1, j1 + 1):
        s = sum(v * v for (mm, q, M), v in out.items() if mm == m)
        worst = max(worst, abs(s - 1.0))
    return worst


print("=" * 74)
print("Clebsch-Gordan coefficients from J^2 diagonalisation")
print("=" * 74)

# ---- validate the construction on a case with an exact known answer ----
print("\n[Validation] 0 (x) 1  ->  J = 1   (expect |CG| = 1 for all q)")
t = cg_from_diagonalisation(0, 1)
for q in (1, 0, -1):
    print(f"   m1=0 q={q:+d} M={q:+d} : CG = {t[(0, q, q)]:+.6f}")
print(f"   max orthonormality error = {check_orthonormality(t, 0, 1):.2e}")

print("\n[Validation] 1 (x) 1  ->  J = 0,1,2")
t11 = cg_from_diagonalisation(1, 1)
print(f"   |J=0, M=0> from (m1=1,q=-1): CG = {t11[(1, -1, 0)]:+.6f}  (expect -0.5774)")
print(f"   max orthonormality error = {check_orthonormality(t11, 1, 1):.2e}")

# ---- the real case: J_l = 2 (lower 3P2), photon 1, J_u = 1 (upper 3S1) ----
print("\n" + "=" * 74)
print("Real case: 6s6p 3P2 (J_l=2)  ->  6s7s 3S1 (J_u=1), photon j=1")
print("=" * 74)
CG = cg_from_diagonalisation(2, 1)
print(f"max orthonormality error = {check_orthonormality(CG, 2, 1):.2e}")

g_u, g_l = 2.0, 1.5   # Lande factors: 3S1 -> 2, 3P2 -> 3/2

rows = []
for m_u in (1, 0, -1):
    for m_l in (2, 1, 0, -1, -2):
        q = m_u - m_l
        if abs(q) > 1:
            continue
        key = (m_l, q, m_u)
        if key not in CG:
            continue
        S = CG[key] ** 2
        delta = g_u * m_u - g_l * m_l
        pol = "pi" if q == 0 else ("sigma+" if q == 1 else "sigma-")
        rows.append((delta, m_u, m_l, q, pol, S))

rows.sort(key=lambda r: -r[0])
tot = sum(r[5] for r in rows)
print(f"\nsum of |CG|^2 over all {len(rows)} allowed transitions = {tot:.6f}")
print(f"(expected 2*J_l+1 = 5 for a complete set)\n")

print(f"{'shift':>7} {'pol':>7} {'m_u':>4} {'m_l':>4} {'q':>3} "
      f"{'|CG|^2':>9} {'rel %':>7}")
print("-" * 56)
for delta, m_u, m_l, q, pol, S in rows:
    print(f"{delta:7.2f} {pol:>7} {m_u:4d} {m_l:4d} {q:3d} {S:9.5f} "
          f"{100 * S / tot:7.2f}")
print("-" * 56)
print(f"{'SUM':>7} {'':>7} {'':>4} {'':>4} {'':>3} {tot:9.5f} {100.0:7.2f}")

# ---- group into the observed pattern ----
print("\nGrouped by polarisation class:")
for name, sel in (("pi  (dm = 0)", lambda r: r[4] == "pi"),
                  ("sigma+ (dm = +1)", lambda r: r[4] == "sigma+"),
                  ("sigma- (dm = -1)", lambda r: r[4] == "sigma-")):
    s = sum(r[5] for r in rows if sel(r))
    print(f"   {name:18s} total |CG|^2 = {s:.5f}  ({100 * s / tot:5.2f} %)")
