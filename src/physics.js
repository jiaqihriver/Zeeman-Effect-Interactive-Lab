/* ==========================================================================
   physics.js — Zeeman effect physics engine (framework-free, testable)
   Hg I 546.1 nm : 6s7s 3S1 (J=1)  ->  6s6p 3P2 (J=2)
   Data: NIST ASD (Hg I), CODATA 2018.
   No dependencies. Works in node (module.exports) and in the browser.
   ========================================================================== */
(function (root, factory) {
  if (typeof module === 'object' && module.exports) module.exports = factory();
  else root.ZeemanPhysics = factory();
})(typeof self !== 'undefined' ? self : this, function () {
  'use strict';

  /* ---------------------------------------------------------------- constants */
  const H_PLANCK = 6.62607015e-34;        // J s      (exact, SI 2019)
  const C_LIGHT  = 2.99792458e8;          // m/s      (exact, SI 2019)
  const K_B      = 1.380649e-23;          // J/K      (exact, SI 2019)
  const MU_B     = 9.2740100783e-24;      // J/T      (CODATA 2018)
  const E_H      = 1.602176634e-19;       // C        (exact, SI 2019)
  const H_C      = H_PLANCK * C_LIGHT;    // J m
  const H_C_CM   = H_PLANCK * C_LIGHT * 100; // J cm  (wavenumber in cm^-1)

  /** 1 Lorentz unit = mu_B * B / (h c) expressed in cm^-1 per tesla. */
  const LORENTZ_PER_T = MU_B / H_C_CM;    // 0.46686... cm^-1 / T

  /* ------------------------------------------------------------ Hg I levels */
  // NIST Atomic Spectra Database, Hg I. Wavenumbers in cm^-1 (vacuum).
  const LEVELS = {
    // ground / lower manifold 6s6p 3P
    P0: { name: '6s6p ³P₀', J: 0, L: 1, S: 1, E: 37644.982, cfg: '6s6p' },
    P1: { name: '6s6p ³P₁', J: 1, L: 1, S: 1, E: 39412.237, cfg: '6s6p' },
    P2: { name: '6s6p ³P₂', J: 2, L: 1, S: 1, E: 44042.909, cfg: '6s6p' },
    // upper 6s7s 3S1  — the 546.1 nm line
    S1: { name: '6s7s ³S₁', J: 1, L: 0, S: 1, E: 62350.325, cfg: '6s7s' }
  };

  const LINE = {
    id: 'HgI-546',
    upper: 'S1',
    lower: 'P2',
    lambda0_nm_air: 546.07,               // standard air value used in labs
    filterLabel: '546.1 nm 干涉滤光片'
  };

  /* --------------------------------------------------------------- Landé g */
  /** Landé g factor: g = 1 + [J(J+1) + S(S+1) - L(L+1)] / (2 J (J+1)) */
  function landeG(L, S, J) {
    if (J === 0) return 0;
    return 1 + (J * (J + 1) + S * (S + 1) - L * (L + 1)) / (2 * J * (J + 1));
  }

  /* ------------------------------------------------------------------------
     Clebsch-Gordan line strengths for the 3S1 -> 3P2 multiplet.

     Values are the squared CG coefficients |⟨J_l m_l ; 1 q | J_u m_u⟩|²,
     obtained independently by diagonalising J² in the uncoupled basis
     (see tools/cg_exact.py).  Each fixed-m_u slice sums to exactly 1, which is
     the orthonormality guarantee for the rotation matrix.

       m_l :  q=-1   q=0    q=+1
       +2  :  0.6    -      -
       +1  :  0.3   0.3     -
        0  :  0.1   0.4    0.1
       -1  :   -    0.3    0.3
       -2  :   -     -     0.6
     ---------------------------------------------------------------------- */
  const CG2 = {
    '2': { '-1': 0.6 },
    '1': { '-1': 0.3, '0': 0.3 },
    '0': { '-1': 0.1, '0': 0.4, '1': 0.1 },
    '-1': { '0': 0.3, '1': 0.3 },
    '-2': { '1': 0.6 }
  };

  function cgStrength(mL, q) {
    const row = CG2[String(mL)];
    if (!row) return 0;
    const v = row[String(q)];
    return v === undefined ? 0 : v;
  }

  /* ------------------------------------------------------------------ spectra */
  const NUD = LINE.upper === 'S1' ? 0 : 0; // vacuum wavenumber of the line
  const NU0 = LEVELS[LINE.upper].E - LEVELS[LINE.lower].E;   // 18307.416 cm^-1
  const LAMBDA0_NM = 1e7 / NU0;                              // 546.227 nm (vacuum)

  const G_UP = landeG(LEVELS.S1.L, LEVELS.S1.S, LEVELS.S1.J);  // 2.0
  const G_LO = landeG(LEVELS.P2.L, LEVELS.P2.S, LEVELS.P2.J);  // 1.5

  /**
   * Build the full set of allowed Zeeman components at field B (tesla).
   * Returns an array of component objects, sorted by shift (descending).
   *
   *  delta = g_u m_u - g_l m_l        (shift in Lorentz units)
   *  dE    = delta * mu_B * B         (J)
   *  nu    = NU0 + dE/(h c)           (cm^-1,  +dE => higher wavenumber)
   *  lambda= 1e7 / nu                 (nm, vacuum)
   */
  function components(B) {
    const out = [];
    const Ju = LEVELS[LINE.upper].J;
    const Jl = LEVELS[LINE.lower].J;
    for (let mU = Ju; mU >= -Ju; mU--) {
      for (let mL = Jl; mL >= -Jl; mL--) {
        const q = mU - mL;                       // Δm  (photon angular momentum)
        if (q < -1 || q > 1) continue;          // electric-dipole selection rule
        const s = cgStrength(mL, q);
        if (s <= 0) continue;
        const delta = G_UP * mU - G_LO * mL;     // in Lorentz units
        const dE = delta * MU_B * B;            // J
        const nu = NU0 + dE / H_C_CM;           // cm^-1
        const lambda = 1e7 / nu;                // nm
        out.push({
          mU, mL, q: mU - mL, delta, strength: s,
          dE_J: dE,
          dE_over_muB: delta * B,               // in units of mu_B*B
          nu_cm: nu,
          dnu_cm: dE / H_C_CM,
          lambda_nm: lambda,
          dlambda_nm: lambda - LAMBDA0_NM,
          pol: q === 0 ? 'pi' : (q > 0 ? 'sigma+' : 'sigma-'),
          // visibility depends on observation direction / polariser, filled later
          visPerp: q === 0 ? 'parallel' : 'perpendicular',
          visPara: q === 0 ? 'dark' : 'circular'
        });
      }
    }
    out.sort((a, b) => b.delta - a.delta);
    return out;
  }

  /**
   * Filter components for a given observation geometry.
   *   view: 'perp'  — observer perpendicular to B
   *         'para'  — observer along B
   *   polariser: 'none' | 'parallel' | 'perpendicular'
   * Returns a new array carrying a `weight` (0..1) for rendering.
   *
   * Physics of polarization:
   *   Δm = 0  (π):  E ∥ B.  Absent when looking along B.
   *   Δm = ±1 (σ):  circularly polarized about B.  When viewed ⊥ to B the
   *                 light is linearly polarised ⊥ to B.  When viewed ALONG B
   *                 it stays circular, and a linear polariser transmits half
   *                 of it irrespective of its orientation — it cannot
   *                 distinguish σ⁺ from σ⁻.
   *
   * A linear polariser at angle φ from B (0 = ∥B, 90 = ⊥B) then transmits,
   * by Malus's law I = I0 cos²φ :
   *   viewed ⊥ B :  π  ->  cos²φ        σ  ->  sin²φ
   *   viewed ∥ B :  π  ->  0            σ  ->  1/2 (independent of φ)
   * The transmission is continuous in φ, so the multiplet fades smoothly as the
   * polariser is rotated, which is what a real bench experiment shows.
   *
   * @param polariser  'none' | 'parallel' | 'perpendicular'  (discrete, kept
   *                   for compatibility)  OR  a number in degrees 0..360.
   *
   * NOTE: never write `opts.polariser || 'none'` to supply a default — 0 is a
   * falsy number and would be swallowed, silently showing the unfiltered
   * multiplet at 0 deg.  Always route through normalisePolariser().
   */

  /** Canonical polariser setting: null (removed) or an angle in degrees 0..360. */
  function normalisePolariser(p) {
    if (p === undefined || p === null) return null;
    if (typeof p === 'string') {
      if (p === 'none') return null;
      return p === 'parallel' ? 0 : 90;
    }
    if (typeof p === 'number' && isFinite(p)) return ((p % 360) + 360) % 360;
    return null;
  }
  function polariserWeight(polariser, isPi) {
    const ang = normalisePolariser(polariser);
    if (ang === null) return 1;
    const phi = ang * Math.PI / 180;
    const c2 = Math.cos(phi) * Math.cos(phi);
    return isPi ? c2 : 1 - c2;          // sin²φ
  }

  /* =============================================== quarter-wave plate (QWP) */
  /**
   * Canonical QWP setting: null (not inserted) or the fast-axis angle in
   * degrees, measured from B in the plane perpendicular to the beam.
   * The plate is a linear retarder: it delays the slow axis by a quarter period,
   * which is exactly lambda/4 — hence "quarter wave".
   */
  function normaliseQWP(q) {
    if (q === undefined || q === null) return null;
    if (typeof q === 'number' && isFinite(q)) return ((q % 180) + 180) % 180;
    return null;
  }

  /**
   * Transmission of one circularly-polarised component through a QWP at
   * fast-axis angle `qw` followed by a linear polariser at angle `pol`.
   *
   * The 2x angle doubling is the whole point of the plate and the reason it is
   * used as a circular-polarisation separator:
   *
   *   Circular light entering a QWP emerges LINEAR, with its polarisation
   *   direction rotated by  TWICE the plate angle:
   *       sigma+  ->  linear at  2*qw
   *       sigma-  ->  linear at  2*qw + 90deg
   *
   * Malus's law then applies to that linear light:
   *       w(sigma+) = cos^2(pol - 2*qw)
   *       w(sigma-) = cos^2(pol - 2*qw - 90deg) = sin^2(pol - 2*qw)
   *
   * Special cases that make the device easy to demonstrate:
   *   qw = pol/2      ->  sigma- fully extinguished (pure sigma+ passes)
   *   qw = pol/2 +45  ->  sigma+ fully extinguished (pure sigma- passes)
   *   qw - pol = 45   ->  both pass at 50 % (the "non-separating" setting)
   *
   * @param q      QWP fast-axis angle in degrees, or null when not inserted
   * @param pol    linear polariser angle in degrees, or null when removed
   * @param qSign  +1 for sigma+, -1 for sigma-  (which circular group)
   * @returns transmission 0..1, or null when the QWP is not in the beam
   */
  function qwpWeight(q, pol, qSign) {
    const qa = normaliseQWP(q);
    if (qa === null) return null;              // plate not inserted
    const pa = normalisePolariser(pol);
    if (pa === null) return 1;                 // no polariser behind: all passes
    const phi = (pa - 2 * qa) * Math.PI / 180;
    const c2 = Math.cos(phi) * Math.cos(phi);
    // the two groups leave the plate 90 deg apart, so they differ by sin^2
    return qSign > 0 ? c2 : 1 - c2;
  }

  /** True when the QWP + polariser fully passes one circular group and kills
   *  the other, i.e. the device is acting as a circular-polarisation separator. */
  function qwpIsSeparator(q, pol) {
    const qa = normaliseQWP(q), pa = normalisePolariser(pol);
    if (qa === null || pa === null) return false;
    const phi = (pa - 2 * qa) * Math.PI / 180;
    const c2 = Math.cos(phi) * Math.cos(phi);
    return c2 < 1e-6 || c2 > 1 - 1e-6;
  }

  function applyGeometry(comps, view, polariser, qwp) {
    const ang = normalisePolariser(polariser);
    const qa = normaliseQWP(qwp);
    const hasPol = ang !== null;
    const hasQ = qa !== null && view === 'para';
    return comps.map(c => {
      const isPi = c.q === 0;
      let w = 1, note = '';
      if (view === 'para') {
        if (isPi) {
          w = 0; note = 'π 沿场方向不辐射';
        } else if (hasQ) {
          /* QWP + polariser: the plate turns the circular light into linear
             light at 2*qw (sigma+) and 2*qw+90 (sigma-); Malus then applies. */
          w = qwpWeight(qa, ang, c.q > 0 ? 1 : -1);
          const turned = ((2 * qa) % 180 + 180) % 180;
          note = `圆偏振 →线偏振 @${turned.toFixed(0)}°（QWP 倍角 2×${qa}°），` +
                 `透射 ${(w * 100).toFixed(0)}%`;
        } else {
          note = '圆偏振';
          // a linear polariser halves circular light, whatever its axis
          if (hasPol) { w = 0.5; note += '，经偏振片减半（无法区分左右旋）'; }
        }
      } else {
        // viewed perpendicular to B: pi is linearly polarised along B,
        // sigma perpendicular to B;  Malus's law applies to the polariser
        w = polariserWeight(ang, isPi);
        if (!hasPol) note = isPi ? '线偏振 ∥B' : '线偏振 ⊥B';
        else if (w < 1e-4) note = '被偏振片滤除';
        else note = isPi ? `透射 ${(w * 100).toFixed(0)}%（∥B 分量）`
                         : `透射 ${(w * 100).toFixed(0)}%（⊥B 分量）`;
      }
      return Object.assign({}, c, {
        weight: w, note,
        polAngle: hasPol ? ((ang % 180) + 180) % 180 : null,
        qwpAngle: hasQ ? qa : null,
        // which of the two circular groups this line belongs to, and which
        // upper sub-level it originates from — used for the group annotation
        group: c.q > 0 ? 'sigma+' : c.q < 0 ? 'sigma-' : 'pi',
        originM: c.mU
      });
    });
  }

  /**
   * Weight below which a component counts as extinguished.
   * Needed because cos^2(phi) reaches 3.7e-33 (not exact 0) at phi = 90 deg.
   */
  const VISIBLE_EPS = 1e-3;

  /** Components that actually transmit, given a geometry. */
  function visibleComponents(comps, view, polariser, qwp) {
    return applyGeometry(comps, view, polariser, qwp)
      .filter(c => c.weight >= VISIBLE_EPS);
  }

  /**
   * Group the sigma components by their circular polarisation sense and report
   * the observable summary used to annotate the spectrum:
   *   { 'sigma+': {n, from, to, lines:[...]}, 'sigma-': {...} }
   * `from`/`to` bracket the group's wavenumber span so the caller can draw a
   * group bracket that never overlaps the neighbouring group.
   */
  function sigmaGroups(comps) {
    const out = {};
    ['sigma+', 'sigma-'].forEach(g => {
      const lines = comps.filter(c => c.group === g)
        .sort((a, b) => a.delta - b.delta);
      out[g] = {
        n: lines.length,
        vis: lines.filter(c => c.weight >= VISIBLE_EPS).length,
        lines,
        from: lines.length ? Math.min.apply(null, lines.map(c => c.delta)) : 0,
        to: lines.length ? Math.max.apply(null, lines.map(c => c.delta)) : 0,
        // which upper sub-levels feed this group
        origins: lines.map(c => c.originM)
      };
    });
    return out;
  }

  /* ------------------------------------------------------- spectral helpers */
  /** Convert a wavenumber to vacuum wavelength in nm. */
  function nuToLambda(nu) { return 1e7 / nu; }
  /** Convert a vacuum wavelength in nm to wavenumber in cm^-1. */
  function lambdaToNu(nm) { return 1e7 / nm; }

  /** Neighbouring-component spacing (mean over adjacent visible components). */
  function meanSpacing(comps) {
    const vis = comps.filter(c => c.weight > 0);
    if (vis.length < 2) return 0;
    let sum = 0, n = 0;
    for (let i = 1; i < vis.length; i++) {
      sum += Math.abs(vis[i].nu_cm - vis[i - 1].nu_cm);
      n++;
    }
    return n ? sum / n : 0;
  }

  /* ============================================================== Zeeman fit */
  /**
   * Ordinary least squares  y = a x + b.
   * Returns slope, intercept, r2, and standard errors.
   */
  function linfit(xs, ys) {
    const n = xs.length;
    if (n < 2) return null;
    let sx = 0, sy = 0, sxx = 0, sxy = 0, syy = 0;
    for (let i = 0; i < n; i++) {
      sx += xs[i]; sy += ys[i];
      sxx += xs[i] * xs[i]; sxy += xs[i] * ys[i]; syy += ys[i] * ys[i];
    }
    const den = n * sxx - sx * sx;
    if (Math.abs(den) < 1e-30) return null;
    const a = (n * sxy - sx * sy) / den;
    const b = (sy - a * sx) / n;
    const ybar = sy / n;
    let ssTot = 0, ssRes = 0;
    for (let i = 0; i < n; i++) {
      const pred = a * xs[i] + b;
      ssRes += (ys[i] - pred) * (ys[i] - pred);
      ssTot += (ys[i] - ybar) * (ys[i] - ybar);
    }
    const r2 = ssTot > 0 ? 1 - ssRes / ssTot : 1;
    // slope standard error
    let se = 0;
    if (n > 2 && ssTot > 0) se = Math.sqrt((ssRes / (n - 2)) / (sxx - sx * sx / n));
    return { slope: a, intercept: b, r2, slopeErr: se, n };
  }

  /* ================================================== magnet: B = f(I) model */
  /**
   * C-core electromagnet with air gap.  B saturates as the core approaches
   * saturation, so B(I) is distinctly non-linear — which is exactly why the
   * experiment requires an independent calibration curve.
   *   B(I) = B_sat * tanh( I / I0 )
   * I0 is chosen so that B(4 A) ~ 0.80 T and B(6 A) ~ 1.05 T.
   */
  const MAGNET = { B_sat: 1.30, I0: 5.20, I_max: 6.0 };
  function fieldFromCurrent(I) {
    return MAGNET.B_sat * Math.tanh(Math.max(0, I) / MAGNET.I0);
  }

  /* ================================================= Fabry–Perot etalon */
  /**
   * Airy transmission of a Fabry–Perot etalon.
   *   T(d) = 1 / (1 + F sin^2(pi d cos(theta) / lambda))
   * F = 4R/(1-R)^2 with R the mirror reflectivity; finesse = pi*sqrt(R)/(1-R).
   *
   * Spacing constraint: the free spectral range FSR = 1/(2 n d) must be much
   * larger than the Zeeman span (4 * 0.4669 * B cm^-1, i.e. 1.87 cm^-1 at 1 T),
   * otherwise successive orders overlap and the multiplet aliases.  With
   * d = 1.0 mm the FSR is 5 cm^-1, a 2.7x margin at full field, and a 4 deg
   * field of view shows ~9 rings.
   */
  const FP = {
    d_mm: 1.0,        // etalon spacing  (FSR = 5 cm^-1)
    n: 1.000277,      // air refractive index at 546 nm
    R: 0.90,          // mirror reflectivity
    f_mm: 250.0,      // lens focal length (imaging)
    get finesse() { return Math.PI * Math.sqrt(this.R) / (1 - this.R); },
    get F() { return 4 * this.R / ((1 - this.R) * (1 - this.R)); },
    /** Free spectral range in cm^-1. */
    get fsr_cm() { return 1 / (2 * this.n * (this.d_mm / 10)); },
    /** Number of bright rings within a field of view thetaMax (radians). */
    ringCount(thetaMax) {
      return this.n * (this.d_mm / 10) * Math.tan(thetaMax) ** 2 / (LAMBDA0_NM * 1e-7);
    }
  };

  /**
   * Airy transmission of a Fabry-Perot etalon at angle theta.
   *
   * The round-trip phase difference is
   *     delta = (2 pi / lambda) * 2 n d cos(theta)
   * and the Airy function is
   *     T = 1 / (1 + F sin^2(delta/2))
   *       = 1 / (1 + F sin^2(pi * 2 n d cos(theta) / lambda)).
   *
   * Note the factor 2: transmission maxima occur when 2 n d cos(theta)/lambda
   * is an INTEGER, which is exactly the resonance condition m lambda = 2 n d
   * cos(theta) used by orderAt() and ringRadius().  Keeping the two consistent
   * is essential — otherwise the rendered rings and the predicted orders
   * disagree (and with a half-integer mAxis they end up on alternating orders).
   */
  function airy(cosTheta, lambdaNm, fp) {
    const d_cm = fp.d_mm / 10;            // mm -> cm
    const lambda_cm = lambdaNm * 1e-7;     // nm -> cm
    const halfPhase = Math.PI * 2 * fp.n * d_cm * cosTheta / lambda_cm;
    const s = Math.sin(halfPhase);
    return 1 / (1 + fp.F * s * s);
  }

  /**
   * Interference order m of a ring observed at angle theta for wavelength lambda.
   * Resonant condition:  m lambda = 2 n d cos(theta)  =>  m = 2 n d cos(theta)/lambda
   * Returns a non-integer "order" whose integer part labels the ring.
   */
  function orderAt(theta, lambdaNm, fp) {
    const d_cm = fp.d_mm / 10;
    const lambda_cm = lambdaNm * 1e-7;
    return 2 * fp.n * d_cm * Math.cos(theta) / lambda_cm;
  }

  /**
   * Ring radius (mm) for order m at wavelength lambda, imaged with focal f.
   * From m lambda = 2 n d cos(theta)  and  r = f tan(theta):
   *    cos(theta) = m lambda / (2 n d)
   *    r = f * tan(arccos(m lambda / 2 n d))
   */
  function ringRadius(m, lambdaNm, fp) {
    const d_cm = fp.d_mm / 10;
    const lambda_cm = lambdaNm * 1e-7;
    const c = m * lambda_cm / (2 * fp.n * d_cm);
    if (c <= 0 || c >= 1) return NaN;
    return fp.f_mm * Math.tan(Math.acos(c));
  }

  /**
   * The quantity a student measures: the squared radius of the OUTERMOST ring
   * (the reddest component, delta = -D/2) as a function of field.
   *
   *   r(m, lambda)^2 = f^2 [ (2 n d nu / m)^2 - 1 ],   nu = 1/lambda
   *
   * Here 2 n d is in cm and nu in cm^-1, so (2 n d nu / m) is dimensionless
   * exactly as required.  No extra scaling factors are needed.
   *
   * The outermost ring has the SMALLEST wavenumber:
   *   nu(B) = nu_0 - (D/2) mu_B B /(h c)      [cm^-1]
   * giving
   *   r^2(B) = f^2 [ (2 n d nu(B) / m)^2 - 1 ]
   * which is quadratic in B; over the experimental range the quadratic term is
   * negligible (1.6e-2 mm^2 out of 1.5e2 mm^2 at 1 T) so a linear fit of r^2
   * versus B is excellent.
   */
  function ringRadiusSqVsB(m, B, D, fp) {
    const d_cm = fp.d_mm / 10;
    const shift = (D / 2) * LORENTZ_PER_T * B;         // cm^-1, red side
    const nu = NU0 - shift;                           // cm^-1
    const x = 2 * fp.n * d_cm * nu / m;               // dimensionless = cos(theta)
    return fp.f_mm * fp.f_mm * (x * x - 1);           // mm^2
  }

  /** Theoretical slope of r^2 vs B (mm^2/T) for the outermost ring. */
  function ringSlope(m, D, fp) {
    const d_cm = fp.d_mm / 10;
    // d(r^2)/dB = 2 f^2 (2nd/m)^2 nu * d(nu)/dB,  d(nu)/dB = -(D/2) LORENTZ_PER_T
    const k = 2 * fp.n * d_cm / m;                    // cm
    const dnu = -(D / 2) * LORENTZ_PER_T;             // cm^-1 per T
    return 2 * fp.f_mm * fp.f_mm * k * k * NU0 * dnu;  // mm^2/T
  }

  /**
   * Invert a measured slope d(r^2)/dB to obtain mu_B.
   *
   *   slope = -2 f^2 K nu_0 (D/2) mu_B /(h c) / 100      [mm^2/T]
   *   K     = (2 n d / m)^2  in cm^-2  ->  x 1e4 to give (cm^-1)^2
   *
   * Solving for mu_B, with all lengths internally in cm:
   */
  function muBFromSlope(slope_mm2_per_T, m, D, fp) {
    const d_cm = fp.d_mm / 10;
    const f_cm = fp.f_mm / 10;
    // slope is in mm^2/T -> cm^2/T  (1 cm^2 = 100 mm^2)
    const s = Math.abs(slope_mm2_per_T) / 100;          // cm^2/T
    const k = 2 * fp.n * d_cm / m;                       // cm
    // |s| = 2 f^2 k^2 nu_0 (D/2) mu_B /(h c)   [cm^2/T]
    const denom = 2 * f_cm * f_cm * k * k * NU0 * (D / 2);
    return s * H_C_CM / denom;
  }

  /** The order D (difference in Lorentz units) for the outermost pair. */
  const OUTER_D = 4;   // delta_+ = +2, delta_- = -2

  /* ======================================================================
     Real-lab measurement protocol
     ----------------------------------------------------------------------
     Step 1 (B = 0): measure the DIAMETERS of the single rings of two
             adjacent interference orders k and k-1.  From the pair
             (D_k, D_{k-1}) one recovers the geometrical constant
             u = lambda0 /(2 n d)  and hence the absolute order k.

     Step 2 (B != 0): the order-(k-1) group now splits into nine rings, one
             per Zeeman component, numbered 1..9 from the INNERMOST ring
             outwards.  Measure the diameters D_p and D_q of a symmetric
             pair (p, q) and extract mu_B.
     ====================================================================== */

  /**
   * Line index 1..9 counted from the INNERMOST ring outwards, and its shift
   * in Lorentz units.
   *
   * A longer wavelength gives a smaller cos(theta) = m lambda/(2nd), hence a
   * larger angle and a larger radius, so the reddest components (delta < 0)
   * form the OUTER rings.  The ordering is therefore
   *
   *   index:1     2     3     4     5     6     7     8     9
   *   delta:-2  -1.5  -1.0  -0.5   0   +0.5  +1.0  +1.5  +2.0
   *
   * i.e.  delta(index) = 2.5 - index  Lorentz units.  Index 5 is the
   * unshifted line (delta = 0) sitting exactly at the zero-field radius, which
   * makes it the natural reference ring.
   */
  const LINE_INDEX_DELTA = [2, 1.5, 1, 0.5, 0, -0.5, -1, -1.5, -2];

  /** Shift in Lorentz units of the line with the given 1..9 index. */
  function deltaOfLineIndex(idx) {
    return LINE_INDEX_DELTA[idx - 1];
  }

  /** Line index 1..9 (inner to outer) of a component with shift delta. */
  function lineIndexOfDelta(delta) {
    const i = LINE_INDEX_DELTA.findIndex(d => Math.abs(d - delta) < 1e-9);
    return i < 0 ? null : i + 1;
  }

  /**
   * The two symmetric line pairs the lab manual uses, given as inner/outer
   * index pairs.  Both are symmetric about delta = 0:
   *   (4, 6): delta = -0.5 / +0.5  ->  the outer pair of pi lines
   *   (3, 7): delta = -1.0 / +1.0  ->  the bright sigma pair
   * `d0` is the magnitude of the shift, in Lorentz units, and the two rings
   * separate by 2*d0 in shift, so the diameter difference grows with d0 and
   * the (3,7) pair is the more sensitive choice.
   */
  const SYMMETRIC_PAIRS = [
    { inner: 4, outer: 6, d0: 0.5, tag: '第 4 / 6 条（π 线对）' },
    { inner: 3, outer: 7, d0: 1.0, tag: '第 3 / 7 条（σ 线对）' }
  ];

  /**
   * Step 1: recover u = lambda0/(2 n d) and the order k from the zero-field
   * diameters of two adjacent orders.
   *
   * With  cos(theta) = m lambda0 /(2 n d) = m u  and  r = f tan(theta):
   *     cos(theta_m) = f / sqrt(f^2 + r_m^2)
   * Adjacent orders differ by exactly one step:
   *     u = cos(theta_k) - cos(theta_{k-1})        (k is the inner ring)
   * and the order then follows algebraically from  k = cos(theta_k) / u.
   *
   * That inversion is exact for perfect readings but far too fragile in
   * practice.  u is a difference of two nearly equal cosines (~2.7e-4), so a
   * ring-position error of only 4e-4 mm shifts u by ~0.2 % and therefore k by
   * SEVERAL WHOLE ORDERS.  Real bench readings always carry that much error, so
   * an algebraic inversion is not usable at the bench.
   *
   * The order is therefore determined by SEARCH: for each candidate order the
   * theoretical pair of radii is compared against the measured pair and the
   * best-fitting order wins.  This uses the measurement the way a numerical fit
   * should, and is stable to a few parts in 10^4.
   *
   * @param Dk       diameter (mm) of the inner ring  (order k)
   * @param Dk1      diameter (mm) of the outer ring  (order k-1)
   * @param lambdaNm central wavelength (nm)
   * @param search   optional {lo, hi} bracket for the candidate order
   * @returns {k, u, uMeasured, cosK, cosK1, residual}
   */
  function orderFromTwoDiameters(Dk, Dk1, lambdaNm, fp, search) {
    fp = fp || FP;
    const f = fp.f_mm;
    const rk = Dk / 2, rk1 = Dk1 / 2;
    if (!(rk > 0 && rk1 > rk)) return null;         // k-1 must be the outer ring

    const cosK = f / Math.sqrt(f * f + rk * rk);
    const cosK1 = f / Math.sqrt(f * f + rk1 * rk1);
    const uMeasured = cosK - cosK1;                // fragile, reported only

    const lo = (search && search.lo) || 2;
    const hi = (search && search.hi) || Math.ceil(2 * fp.n * (fp.d_mm / 10) * 1e7 / lambdaNm);
    let best = null;
    for (let k = Math.max(2, lo); k <= hi; k++) {
      const tk = ringRadius(k, lambdaNm, fp);
      const tk1 = ringRadius(k - 1, lambdaNm, fp);
      if (!isFinite(tk) || !isFinite(tk1)) continue;
      // relative misfit on BOTH radii, so the larger ring is not favoured
      const res = Math.hypot((rk - tk) / tk, (rk1 - tk1) / tk1);
      if (!best || res < best.res) best = { k, res, tk, tk1 };
    }
    if (!best) return null;

    const u = lambdaNm * 1e-7 / (2 * fp.n * (fp.d_mm / 10));
    return {
      k: best.k,
      kExact: uMeasured > 0 ? cosK / uMeasured : best.k,
      u,
      uMeasured,
      uRelError: uMeasured > 0 ? Math.abs(uMeasured - u) / u : NaN,
      cosK, cosK1,
      residual: best.res
    };
  }

  /**
   * Step 2: mu_B from the diameters of a symmetric line pair at one order.
   *
   * For a pair at shifts +-d0 (Lorentz units) the wavenumbers are
   *     nu_p = nu0 + d0 L,   nu_q = nu0 - d0 L,   L = mu_B B /(h c)
   * and since  r^2 = f^2 ( K^2 nu^2 - 1 )  with  K = 2 n d / m,
   *     r_p^2 - r_q^2 = f^2 K^2 [ (nu0+d0L)^2 - (nu0-d0L)^2 ]
   *                  = 4 f^2 K^2 nu0 d0 L
   * With DIAMETERS D = 2r the same relation carries a factor 4:
   *     D_p^2 - D_q^2 = 16 f^2 K^2 nu0 d0 mu_B B /(h c)
   * hence
   *     mu_B = (D_p^2 - D_q^2) h c / (16 f^2 K^2 nu0 d0 B)
   *
   * Linear in B to first order, so a plot of (D_p^2 - D_q^2) against B is a
   * straight line through the origin whose slope gives mu_B.
   *
   * @param Dp     diameter (mm) of the OUTER ring of the pair (larger radius:
   *               the redder, negative-delta component)
   * @param Dq     diameter (mm) of the INNER ring (positive delta)
   * @param d0     shift magnitude in Lorentz units (0.5 for lines 4/6,
   *               1.0 for lines 3/7)
   * @param m      interference order of the measured group
   */
  function muBFromDiameterPair(Dp, Dq, d0, m, B, fp) {
    fp = fp || FP;
    if (!(B > 0) || !(m > 0) || !(d0 > 0)) return NaN;
    const d_cm = fp.d_mm / 10;
    const f_cm = fp.f_mm / 10;
    const K = 2 * fp.n * d_cm / m;                 // cm
    // diameters are mm, so D^2 -> cm^2 needs /100
    const dSq = (Dp * Dp - Dq * Dq) / 100;         // cm^2
    return dSq * H_C_CM / (16 * f_cm * f_cm * K * K * NU0 * d0 * B);
  }

  /** Slope of (D_p^2 - D_q^2) against B, in cm^2/T, for a symmetric pair. */
  function diameterPairSlopeScale(m, d0, fp) {
    fp = fp || FP;
    const d_cm = fp.d_mm / 10;
    const f_cm = fp.f_mm / 10;
    const K = 2 * fp.n * d_cm / m;
    // slope = 16 f^2 K^2 nu0 d0 mu_B /(h c)  ->  invert for mu_B
    return 16 * f_cm * f_cm * K * K * NU0 * d0;
  }

  /**
   * Predicted DIAMETER (mm) of a split line, given the field, the order, and
   * the 1..9 line index counted from the innermost ring outwards.
   */
  function splitDiameter(m, lineIndex, B, fp) {
    fp = fp || FP;
    const delta = deltaOfLineIndex(lineIndex);
    if (delta === undefined) return NaN;
    const nu = NU0 - delta * LORENTZ_PER_T * B;   // cm^-1
    const r = ringRadius(m, nuToLambda(nu), fp);
    return 2 * r;
  }

  /**
   * Predicted DIAMETER (mm) of the un-split zero-field ring of order m.
   * Independent of B, because all nine components coincide when B = 0.
   */
  function zeroFieldDiameter(m, fp) {
    fp = fp || FP;
    return 2 * ringRadius(m, LAMBDA0_NM, fp);
  }

  /**
   * Highest interference order visible within a field of view thetaMax.
   *
   * The resonance condition is  m = 2 n d cos(theta) / lambda.  The LARGEST
   * order mAxis = 2 n d / lambda sits at theta = 0 (the very centre); orders
   * below it appear at progressively larger angles, with
   *   r_k = f tan(theta_k),  cos(theta_k) = m_k / mAxis.
   * A ring of order m is therefore inside the field when
   *   m >= mAxis cos(thetaMax).
   * This returns that LOWER bound: orders from minVisibleOrder up to
   * floor(mAxis) are all visible (minVisibleOrder = the outermost ring that
   * still fits inside the field of view).
   */
  function minVisibleOrder(thetaMax, fp) {
    fp = fp || FP;
    const d_cm = fp.d_mm / 10;
    const lambda_cm = LAMBDA0_NM * 1e-7;
    const mAxis = 2 * fp.n * d_cm / lambda_cm;
    return Math.ceil(mAxis * Math.cos(thetaMax));
  }

  /** Highest order present at all (sits at the centre, theta = 0). */
  function maxOrder(fp) {
    fp = fp || FP;
    const d_cm = fp.d_mm / 10;
    const lambda_cm = LAMBDA0_NM * 1e-7;
    return Math.floor(2 * fp.n * d_cm / lambda_cm);
  }

  /** Order of the ring at ring-index k counted outward from the centre. */
  function orderAtIndex(k, thetaMax, fp) {
    fp = fp || FP;
    return Math.max(minVisibleOrder(thetaMax, fp), maxOrder(fp) - k);
  }

  /** A sensible default order to measure: a few rings in from the axis. */
  function suggestOrder(thetaMax, fp) {
    fp = fp || FP;
    return orderAtIndex(4, thetaMax, fp);
  }

  /* ------------------------------------------------------------------ export */
  return {
    // constants
    H_PLANCK, C_LIGHT, K_B, MU_B, E_H, H_C, H_C_CM, LORENTZ_PER_T,
    // atomic data
    LEVELS, LINE, NU0, LAMBDA0_NM, G_UP, G_LO, landeG, CG2, cgStrength,
    // model
    components, applyGeometry, polariserWeight, normalisePolariser,
    normaliseQWP, qwpWeight, qwpIsSeparator, sigmaGroups,
    visibleComponents, VISIBLE_EPS, nuToLambda, lambdaToNu, meanSpacing,
    linfit, MAGNET, fieldFromCurrent,
    FP, airy, orderAt, ringRadius, ringRadiusSqVsB, ringSlope,
    muBFromSlope, OUTER_D,
    LINE_INDEX_DELTA, deltaOfLineIndex, lineIndexOfDelta, SYMMETRIC_PAIRS,
    orderFromTwoDiameters, muBFromDiameterPair, diameterPairSlopeScale,
    splitDiameter, zeroFieldDiameter,
    maxOrder, minVisibleOrder, orderAtIndex, suggestOrder
  };
});
