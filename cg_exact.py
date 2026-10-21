"""
Authoritative Clebsch-Gordan coefficients for the Hg 546.1 nm Zeeman multiplet,
obtained by diagonalising J^2 = J1^2 + J2^2 + 2 J1.J2 in the uncoupled basis
|j1 m1, j2 m2>.

Emission: upper 6s7s 3S1 (J=1)  ->  lower 6s6p 3P2 (J=2)
The photon carries j2 = 1, q = m_u - m_l.

WHY THIS IS DONE PER FIXED-M SUBSPACE
-------------------------------------
J1^2, J2^2 and J1zJ2z are diagonal in the uncoupled basis, so J_z is diagonal
too.  The eigenvalues J(J+1) are therefore DISTINCT within a fixed-M subspace
(each allowed J contributes exactly one state to that M).  Diagonalising each
fixed-M block separately makes the eigenbasis unique.

The earlier revision of this file diagonalised the full matrix and projected
onto definite M with `block @ (block.T @ amp)`.  That is ill-defined: inside a
degenerate J block `numpy.linalg.eigh` returns an ARBITRARY orthonormal basis,
so the projection returned basis-dependent numbers.  Its own self-test caught
this -- it reported |<1,1;1,-1|0,0>| = +0.408 instead of -0.577 and an
orthonormality error of 2.00 -- but the printed diagnostics were easy to
overlook while the misleading coefficients sat in `physics.js` comments.

The checks in main() are therefore hard failures (assert), not printouts: if
this file ever disagrees with `CG2` in src/physics.js, it must exit non-zero.

NOTE: Python strips `assert` statements under `python -O`, which would turn
every one of those checks into a no-op.  main() therefore refuses to run
under -O rather than silently passing.
"""
import sys

import numpy as np

np.set_printoptions(suppress=True)

TOL = 1e-9


def fixed_m_subspace(j1, j2, M):
    """Uncoupled basis states |j1 m1, j2 m2> with m1 + m2 == M."""
    return [(m1, m2) for m1 in range(-j1, j1 + 1)
            for m2 in range(-j2, j2 + 1) if m1 + m2 == M]


def build_J2(j1, j2, M):
    """Matrix of J^2 inside one fixed-M subspace, in that subspace's basis."""
    sub = fixed_m_subspace(j1, j2, M)
    idx = {k: a for a, k in enumerate(sub)}
    n = len(sub)
    H = np.zeros((n, n))
    for (m1, m2), a in idx.items():
        # J1^2 + J2^2 + 2 J1z J2z
        H[a, a] = j1 * (j1 + 1) + j2 * (j2 + 1) + 2.0 * m1 * m2
        # + J1+J2- + J1-J2+   (no factor 1/2: the 2 in 2 J1.J2 already cancels it)
        b = idx.get((m1 + 1, m2 - 1))
        if b is not None:
            H[a, b] += (np.sqrt((j1 - m1) * (j1 + m1 + 1))
                        * np.sqrt((j2 + m2) * (j2 - m2 + 1)))
        b = idx.get((m1 - 1, m2 + 1))
        if b is not None:
            H[a, b] += (np.sqrt((j1 + m1) * (j1 - m1 + 1))
                        * np.sqrt((j2 - m2) * (j2 + m2 + 1)))
    return H, sub, idx


def cg_from_diagonalisation(j1, j2):
    """
    Return dict[(J, m1, q, M)] = <j1 m1; j2 q | J M> for J = |j1-j2| ... j1+j2.

    J is part of the key on purpose.  Every J contributes exactly one state to
    each fixed-M subspace, so those states share the same (m1, m2) labels; a key
    of (m1, q, M) alone makes the different J overwrite one another and the
    lower-J coefficients are silently lost.

    Sign convention is fixed per-J by making the largest-magnitude component
    positive; only |CG|^2 is physically meaningful and is what physics.js uses.
    """
    out = {}
    for J in range(abs(j1 - j2), j1 + j2 + 1):
        for M in range(-J, J + 1):
            H, sub, idx = build_J2(j1, j2, M)
            H = 0.5 * (H + H.T)          # kill round-off asymmetry
            ev, evec = np.linalg.eigh(H)

            # the block must contain J(J+1) exactly once
            hit = np.where(np.abs(ev - J * (J + 1)) < 1e-8)[0]
            assert len(hit) == 1, (
                f"J={J} M={M}: expected exactly one eigenvalue {J*(J+1)}, got {ev}")

            v = evec[:, hit[0]]
            v = v / np.linalg.norm(v)
            if v[np.argmax(np.abs(v))] < 0:
                v = -v
            for (m1, m2), a in idx.items():
                if abs(v[a]) > 1e-9:
                    out[(J, m1, m2, M)] = v[a]
    return out


def check_orthonormality(out, j1, j2, J):
    """
    Column-sum check at fixed J and fixed M = m1 + q, summing over q.

    For fixed J == j2 and fixed M, every lower-level m1 that can reach it
    contributes exactly one q = M - m1, so the sum of |CG|^2 over those m1 is
    the norm of one column of the CG matrix and must be 1.

    This is the invariant physics.js relies on ("each fixed-m_u slice sums to
    exactly 1").  Summing over m_l at fixed m_u is equivalent, since m1 and M
    determine q; summing over m_l at fixed m1 instead is NOT (for 2x1 at J=1
    that gives 0.6).
    """
    if J != j2:
        raise ValueError("column-sum invariant is defined for J == j2 only")
    worst = 0.0
    for M in range(-J, J + 1):
        s = sum(v * v for (JJ, mm, q, MM), v in out.items()
                if JJ == J and MM == M)
        worst = max(worst, abs(s - 1.0))
    return worst


def main():
    if not __debug__:
        # `python -O` strips asserts, which would turn every check below into
        # a silent no-op.  Fail loudly instead of reporting a false pass.
        raise SystemExit(
            "refusing to run under python -O: the assert-based checks in this "
            "script would be stripped, so the validation would be meaningless")

    print("=" * 74)
    print("Clebsch-Gordan coefficients from J^2 diagonalisation (per fixed M)")
    print("=" * 74)

    # ---- gate 1: the construction must reproduce textbook 1 x 1 values ----
    # These are the numbers the previous (broken) revision got wrong.
    print("\n[Validation] 1 x 1 coupling vs closed-form values")
    t11 = cg_from_diagonalisation(1, 1)
    expect = {
        (0, 1, -1, 0): -1 / np.sqrt(3),     # singlet
        (0, 0, 0, 0): +1 / np.sqrt(3),
        (0, -1, 1, 0): -1 / np.sqrt(3),
        (1, 1, 0, 1): -1 / np.sqrt(2),      # triplet
        (1, 0, 1, 1): +1 / np.sqrt(2),
        (2, 1, -1, 0): 1 / np.sqrt(6),      # highest triplet
        (2, 0, 0, 0): np.sqrt(2 / 3),
    }
    for key, want in expect.items():
        got = t11.get(key)
        assert got is not None, f"missing coefficient {key}"
        assert abs(abs(got) - abs(want)) < TOL, \
            f"{key}: |CG| = {abs(got):.8f}, expected {abs(want):.8f}"
        print(f"   <{key[0]},{key[1]};1,{key[2]}|{key[0]},{key[3]}> "
              f"= {got:+.6f}   OK")
    e = check_orthonormality(t11, 1, 1, 1)
    assert e < TOL, f"1x1 row-sum error {e:.3e}"
    print("   column-sum invariant OK at J = j2 = 1")

    # ---- gate 2: 0 x 1 must be trivial ----
    t01 = cg_from_diagonalisation(0, 1)
    for q in (1, 0, -1):
        assert abs(abs(t01[(1, 0, q, q)]) - 1.0) < TOL, f"0x1 q={q} is not 1"
    print("\n[Validation] 0 x 1 -> J = 1 : all |CG| = 1   OK")

    # ---- the real case: J_l = 2 (lower 3P2), photon 1, J_u = 1 (upper 3S1) --
    print("\n" + "=" * 74)
    print("Real case: 6s6p 3P2 (J_l=2)  ->  6s7s 3S1 (J_u=1), photon j=1")
    print("=" * 74)
    CG = cg_from_diagonalisation(2, 1)
    worst = check_orthonormality(CG, 2, 1, 1)      # J_u = 1 only
    assert worst < TOL, f"orthonormality error {worst:.3e}"
    print(f"max orthonormality error = {worst:.2e}")

    g_u, g_l = 2.0, 1.5   # Lande factors: 3S1 -> 2, 3P2 -> 3/2
    J_U = 1

    rows = []
    for m_u in (1, 0, -1):
        for m_l in (2, 1, 0, -1, -2):
            q = m_u - m_l
            if abs(q) > 1:
                continue
            key = (J_U, m_l, q, m_u)
            if key not in CG:
                continue
            S = CG[key] ** 2
            delta = g_u * m_u - g_l * m_l
            pol = "pi" if q == 0 else ("sigma+" if q == 1 else "sigma-")
            rows.append((delta, m_u, m_l, q, pol, S))

    rows.sort(key=lambda r: -r[0])
    tot = sum(r[5] for r in rows)
    print(f"\nsum of |CG|^2 over all {len(rows)} allowed transitions = {tot:.6f}")
    print("(expected 2*J_u+1 = 3 for a complete multiplet)\n")
    assert len(rows) == 9, f"expected 9 allowed transitions, got {len(rows)}"
    assert abs(tot - 3.0) < TOL, f"sum |CG|^2 = {tot}, expected 3"

    print(f"{'shift':>7} {'pol':>7} {'m_u':>4} {'m_l':>4} {'q':>3} "
          f"{'|CG|^2':>9} {'rel %':>7}")
    print("-" * 56)
    for delta, m_u, m_l, q, pol, S in rows:
        print(f"{delta:7.2f} {pol:>7} {m_u:4d} {m_l:4d} {q:3d} {S:9.5f} "
              f"{100 * S / tot:7.2f}")
    print("-" * 56)
    print(f"{'SUM':>7} {'':>7} {'':>4} {'':>4} {'':>3} {tot:9.5f} {100.0:7.2f}")

    print("\nGrouped by polarisation class:")
    for name, sel in (("pi  (dm = 0)", lambda r: r[4] == "pi"),
                      ("sigma+ (dm = +1)", lambda r: r[4] == "sigma+"),
                      ("sigma- (dm = -1)", lambda r: r[4] == "sigma-")):
        s = sum(r[5] for r in rows if sel(r))
        print(f"   {name:18s} total |CG|^2 = {s:.5f}  ({100 * s / tot:5.2f} %)")

    # ---- gate 3: must agree with the CG2 table hard-coded in physics.js ----
    # Read the table by LOADING physics.js, not by pattern-matching its source:
    # a regex over JS literal syntax is brittle and would itself be a source of
    # false confidence.  physics.js is a UMD module, so node can require() it.
    print("\n[Cross-check] CG2 table in src/physics.js")
    import os
    import subprocess
    import json

    here = os.path.dirname(os.path.abspath(__file__))
    js_path = os.path.join(here, "src", "physics.js")
    node = os.environ.get("NODE", "node")
    probe = (
        "const P=require(process.argv[1]);"
        "const o={};"
        "for(const m of Object.keys(P.CG2)){"
        "  o[m]={};"
        "  for(const q of Object.keys(P.CG2[m])) o[m][q]=P.cgStrength(Number(m),Number(q));"
        "}"
        "process.stdout.write(JSON.stringify(o));"
    )
    try:
        res = subprocess.run([node, "-e", probe, js_path],
                             capture_output=True, text=True, timeout=60)
    except (OSError, subprocess.SubprocessError) as exc:
        raise SystemExit(f"could not run node to load physics.js: {exc}")
    if res.returncode != 0:
        raise SystemExit(f"node failed to load physics.js:\n{res.stderr}")
    cg2 = json.loads(res.stdout)

    mismatches = []
    for delta, m_u, m_l, q, pol, S in rows:
        js_val = cg2.get(str(m_l), {}).get(str(q))
        if js_val is None or abs(js_val - S) > TOL:
            mismatches.append((m_l, q, S, js_val))
    if mismatches:
        for m_l, q, exact, js_val in mismatches:
            print(f"   MISMATCH m_l={m_l:+d} q={q:+d}: exact {exact:.6f} "
                  f"vs CG2 {js_val}")
        raise SystemExit("CG2 in physics.js disagrees with the exact result")
    print("   all 9 entries agree with CG2   OK")


if __name__ == "__main__":
    main()