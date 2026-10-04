/* ==========================================================================
   measure.js — the real bench protocol for extracting the Bohr magneton.

   The laboratory procedure (and the one modelled here) is a TWO-STEP
   measurement.  It is deliberately NOT "track one ring versus B", because
   the etalon geometry contributes a large B-independent term that has to be
   removed before mu_B can be isolated.

   Step 1 — fix the geometry, with B = 0
       All nine Zeeman components coincide, so each interference order shows a
       single ring.  Measure the DIAMETERS of two adjacent orders:
           D_k  (the inner ring)   and   D_{k-1}  (the next ring out)
       From the pair one recovers  u = lambda0 /(2 n d)  and hence the absolute
       interference order k.  This is what fixes the 1/m^2 factor in step 2.
       Equivalent variant: keep the field on and use line 5 of each order.
       Line 5 has delta = 0, so it is field-independent and gives exactly the
       same two reference diameters (verified to 0.00e+0 mm).

   Step 2 — measure a symmetric split pair, with B != 0
       Order (k-1) now shows nine rings, numbered 1..9 from the INNERMOST
       outwards (longer wavelength -> larger radius -> further out, so the
       reddest component is ring 1).  Measure the diameters of a pair that is
       symmetric about delta = 0:
           lines 4 and 6  (delta = -0.5 / +0.5, the outer pi pair),  or
           lines 3 and 7  (delta = -1.0 / +1.0, the bright sigma pair)
       Because
           D_p^2 - D_q^2 = 16 f^2 K^2 nu0 d0 mu_B B /(h c),   K = 2 n d / m,
       the diameter-squared difference is linear in B through the origin, and
       its slope gives mu_B.

   Depends on physics.js.
   ========================================================================== */
(function (root) {
  'use strict';
  const P = root.ZeemanPhysics;
  const R = root.ZeemanRender;

  /** Deterministic PRNG so a given seed always yields the same data set. */
  function mulberry32(seed) {
    let a = seed >>> 0;
    return function () {
      a |= 0; a = (a + 0x6D2B79F5) | 0;
      let t = Math.imul(a ^ (a >>> 15), 1 | a);
      t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
      return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
    };
  }

  /** Box–Muller normal deviate. */
  function gauss(rnd) {
    let u = 0, v = 0;
    while (u === 0) u = rnd();
    while (v === 0) v = rnd();
    return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * v);
  }

  /** State of one measurement session. */
  function Session() {
    this.reset();
    /* The bench observes PERPENDICULAR to B with the polariser removed, so all
       nine components are visible and the line numbering 1..9 is meaningful. */
    this.geom = { view: 'perp', polariser: null, qwp: null };
  }

  Session.prototype.reset = function () {
    this.ref = null;              // step-1 result
    this.rows = [];               // step-2 points
    /* Default to the 3/7 (sigma) pair: it is the brighter pair, so it stays
       resolvable down to lower fields, and it is what the UI highlights. */
    this.pairIdx = 1;
    this.noiseMm = 0.002;         // 1-sigma diameter reading error, mm
    this.seed = 20261004;
    this.rnd = mulberry32(this.seed);
  };

  /** The symmetric pair currently selected. */
  Session.prototype.pair = function () { return P.SYMMETRIC_PAIRS[this.pairIdx]; };
  Session.prototype.setPair = function (i) {
    this.pairIdx = Math.max(0, Math.min(P.SYMMETRIC_PAIRS.length - 1, i | 0));
    return this.pair();
  };

  /** The order whose split rings are measured in step 2 (k-1 from step 1). */
  Session.prototype.splitOrder = function () {
    if (!this.ref) return null;
    return this.ref.k - 1;
  };

  /** Apply the configured reading error to a diameter. */
  Session.prototype.jitter = function (D) {
    if (!(this.noiseMm > 0)) return D;
    return D + this.noiseMm * gauss(this.rnd);
  };

  /* ------------------------------------------------------------------ step 1 */

  /**
   * Step 1: the two zero-field reference rings.
   *
   * @param Dk  diameter (mm) of the INNER ring, interference order k
   * @param Dk1 diameter (mm) of the next ring out, order k-1
   * @returns the order solution, or {ok:false, reason} if the pair is
   *          geometrically impossible (e.g. D_{k-1} <= D_k)
   */
  Session.prototype.setReference = function (Dk, Dk1, source, search) {
    const sol = P.orderFromTwoDiameters(Dk, Dk1, P.LAMBDA0_NM, P.FP, search);
    if (!sol) {
      this.ref = null;
      return { ok: false, reason: '两个环直径不构成相邻级次：需满足 D_{k−1} > D_k。' };
    }
    // the search should land on an order whose pair fits the field of view
    if (!(sol.k >= 2) || sol.residual > 5e-3) {
      this.ref = null;
      return {
        ok: false,
        reason: '在视场范围内找不到与这两个环直径吻合的相邻级次' +
                '（最佳拟合优度 ' + (sol.residual * 100).toFixed(2) + '%）。请检查视场或选更大的 k。'
      };
    }
    this.ref = {
      Dk, Dk1, u: sol.u, uMeasured: sol.uMeasured, uRelError: sol.uRelError,
      k: sol.k, kExact: sol.kExact,
      cosK: sol.cosK, cosK1: sol.cosK1, residual: sol.residual,
      source: source || 'B = 0'
    };
    return { ok: true, ref: this.ref };
  };

  /**
   * Measure the two reference rings straight from the rendered pattern at
   * B = 0, where every order shows exactly one ring.  Rings are located by
   * matching the detected peaks against the theoretical positions, so this
   * exercises the same code path a student would use by eye.
   */
  Session.prototype.measureReference = function (orderGuess, opts) {
    opts = opts || {};
    const fp = P.FP;
    const thetaMax = (opts.thetaMaxDeg || 4) * Math.PI / 180;
    const prof = R.radialProfile(0, Object.assign(
      { N: opts.N || 12000, thetaMaxDeg: opts.thetaMaxDeg || 4 }, this.geom));
    const rings = R.findRings(prof, { minProm: opts.minProm || 0.015 });
    if (rings.length < 2) {
      return { ok: false, reason: '未识别出足够多的干涉环（B = 0 时每级只有一个环）。请增大视场。' };
    }
    const k = Math.round(orderGuess);
    const rk = P.zeroFieldDiameter(k, fp) / 2;
    const rk1 = P.zeroFieldDiameter(k - 1, fp) / 2;
    if (!isFinite(rk) || !isFinite(rk1)) {
      return { ok: false, reason: '级次 ' + k + ' 不在当前视场内，请调整视场或选更小的 k。' };
    }
    const a = nearestRing(rings, rk), b = nearestRing(rings, rk1);
    if (!a || !b) {
      return { ok: false, reason: '未能定位 k 或 k−1 级的环。' };
    }
    // the two rings must be distinct detections
    if (a === b) {
      return { ok: false, reason: 'k 与 k−1 级被识别成同一个环，请检查视场。' };
    }
    return this.setReference(this.jitter(2 * a.r), this.jitter(2 * b.r), 'B = 0',
      { lo: P.minVisibleOrder(thetaMax, fp), hi: P.maxOrder(fp) });
  };

  /* ------------------------------------------------------------------ step 2 */

  /**
   * Step 2: measure one symmetric pair at the current field.
   *
   * The nine rings of the selected order are located by matching the detected
   * profile against the theoretical radii, then the two symmetric ones are
   * read off.  A match is only accepted if it lands closer to its own
   * prediction than to any neighbouring ring, which is what stops a merged or
   * aliased ring from being silently accepted.
   *
   * @param B   field in tesla
   * @param I   coil current in amperes (recorded only)
   */
  Session.prototype.measurePair = function (B, I, opts) {
    opts = opts || {};
    const fp = P.FP;
    if (!this.ref) {
      return { ok: false, reason: '请先完成第 1 步：测量 B = 0 时 k 级与 k−1 级的环直径。' };
    }
    if (!(B > 1e-4)) {
      return { ok: false, reason: '请先加磁场再测量分裂环（B = 0 时 9 条谱线重合）。' };
    }
    const m = this.splitOrder();
    const pr = this.pair();

    // predicted radii of all nine lines of this order
    const pred = [];
    for (let i = 1; i <= 9; i++) {
      pred.push({ idx: i, r: P.splitDiameter(m, i, B, fp) / 2 });
    }
    if (pred.some(p => !isFinite(p.r))) {
      return { ok: false, reason: '级次 m = ' + m + ' 在当前场下的部分谱线无实数半径，请选更大的 k 或增大视场。' };
    }
    // minimum separation anywhere in the pattern: sets the match tolerance
    const sorted = pred.map(p => p.r).sort((a, b) => a - b);
    let minSep = Infinity;
    for (let i = 1; i < sorted.length; i++) minSep = Math.min(minSep, sorted[i] - sorted[i - 1]);
    const tol = 0.45 * minSep;

    // Measure the two rings directly from the profile.  Only a narrow window
    // around each predicted radius is searched, which is both faster and
    // stricter than scanning every ring: a peak belonging to a neighbouring
    // line or a neighbouring interference order can never be picked up.
    const N = opts.N || 20000;
    const rMaxMm = fp.f_mm * Math.tan((opts.thetaMaxDeg || 4) * Math.PI / 180);
    const tolMm = 0.45 * minSep;
    const grab = (targetR) => {
      const prof = R.radialProfile(B, Object.assign(
        { N, thetaMaxDeg: opts.thetaMaxDeg || 4 }, this.geom));
      const idx = targetR / rMaxMm * (prof.N - 1);
      const half = Math.max(3, Math.ceil(tolMm / rMaxMm * (prof.N - 1)));
      const lo = Math.max(1, Math.floor(idx - half));
      const hi = Math.min(prof.N - 2, Math.ceil(idx + half));
      const { prof: v } = prof;
      let best = -1, bi = -1;
      for (let i = lo; i <= hi; i++) {
        if (v[i] > v[i - 1] && v[i] >= v[i + 1] && v[i] > best) { best = v[i]; bi = i; }
      }
      if (bi < 0) return null;
      // parabolic refinement for sub-sample precision
      const a = v[bi - 1], c = v[bi + 1], den = a - 2 * best + c;
      const shift = Math.abs(den) > 1e-18 ? 0.5 * (a - c) / den : 0;
      const r = (bi + Math.max(-1, Math.min(1, shift))) / (prof.N - 1) * rMaxMm;
      return { r, I: best };
    };

    const inner = grab(pred[pr.inner - 1].r);
    const outer = grab(pred[pr.outer - 1].r);
    const miss = [];
    if (!inner) miss.push(pr.inner);
    if (!outer) miss.push(pr.outer);
    if (miss.length) {
      return {
        ok: false,
        reason: '第 ' + miss.join(' 与 ') + ' 条谱线的环未能定位（预测 r = ' +
                pred[miss[0] - 1].r.toFixed(4) + ' mm）。该组相邻线间距仅 ' +
                (minSep * 1000).toFixed(1) + ' μm，低于当前分辨率——请提高 B。'
      };
    }
    inner.predicted = pred[pr.inner - 1].r;
    outer.predicted = pred[pr.outer - 1].r;

    const Di = this.jitter(2 * inner.r);     // inner ring, positive delta
    const Do = this.jitter(2 * outer.r);     // outer ring, negative delta
    const dSq = Do * Do - Di * Di;           // mm^2, positive by construction
    const muB = P.muBFromDiameterPair(Do, Di, pr.d0, m, B, fp);

    const row = {
      B, I, m, pairIdx: this.pairIdx, d0: pr.d0,
      Douter: Do, Dinner: Di,
      rOuter: outer.r, rInner: inner.r,
      DSq: dSq,
      muB,
      driftInner: Math.abs(inner.r - pred[pr.inner - 1].r),
      driftOuter: Math.abs(outer.r - pred[pr.outer - 1].r),
      sepMeasured: outer.r - inner.r,
      sepTheory: pred[pr.outer - 1].r - pred[pr.inner - 1].r
    };
    this.rows.push(row);
    return { ok: true, row };
  };

  /**
   * Run a whole sweep of fields for the selected pair, always including a
   * B = 0 reference first.  Returns a per-point log so the UI can report
   * where the pattern became unresolvable.
   */
  Session.prototype.sweep = function (Bvalues, orderGuess, opts) {
    opts = opts || {};
    this.rows = [];
    const log = [];
    const r0 = this.measureReference(orderGuess, opts);
    log.push({ step: 1, ...r0 });
    if (!r0.ok) return { ok: false, log };
    Bvalues.forEach(B => {
      const r = this.measurePair(B, 0, opts);
      log.push({ step: 2, B, ...r });
    });
    return { ok: true, log };
  };

  /* --------------------------------------------------------------- analysis */

  /**
   * Fit (D_p^2 - D_q^2) against B.  The relation is exactly linear through the
   * origin, so the fit is a single-parameter slope:
   *     slope = sum(B y) / sum(B^2)
   * with slope in cm^2/T, and then
   *     mu_B = slope * h c / (16 f^2 K^2 nu0 d0)
   */
  Session.prototype.analyse = function () {
    if (!this.ref || this.rows.length < 2) return null;
    const m = this.splitOrder();
    const pr = this.pair();
    let sxy = 0, sxx = 0, sy = 0, n = 0;
    const xs = [], ys = [];
    this.rows.forEach(r => {
      if (!(r.B > 0)) return;
      const y = r.DSq / 100;          // mm^2 -> cm^2
      xs.push(r.B); ys.push(y);
      sxy += r.B * y; sxx += r.B * r.B; sy += y; n++;
    });
    if (n < 2) return null;
    const slope_cm2 = sxy / sxx;                       // cm^2/T
    const scale = P.diameterPairSlopeScale(m, pr.d0, P.FP);
    const muB = slope_cm2 * P.H_C_CM / scale;
    // through-origin residuals -> R^2 against the unconstrained best fit
    const interceptFree = P.linfit(xs, ys);
    let ssRes = 0, ssTot = 0;
    const ybar = sy / n;
    this.rows.forEach(r => {
      if (!(r.B > 0)) return;
      const y = r.DSq / 100;
      ssRes += Math.pow(y - slope_cm2 * r.B, 2);
      ssTot += Math.pow(y - ybar, 2);
    });
    const r2 = ssTot > 0 ? 1 - ssRes / ssTot : 0;
    // standard error of the slope for the through-origin fit
    const sigma2 = ssRes / (n - 1);
    const slopeErr = Math.sqrt(sigma2 / sxx);
    const muBErr = Math.abs(muB * slopeErr / slope_cm2);
    const perPoint = this.rows.filter(r => r.B > 0).map(r => r.muB);
    const mean = perPoint.reduce((a, b) => a + b, 0) / perPoint.length;
    return {
      fit: { slope: slope_cm2, slopeErr, intercept: 0, r2, slopeFree: interceptFree.slope },
      m, pair: pr, muB, muBErr,
      meanPerPoint: mean,
      spreadPerPoint: perPoint.length > 1
        ? Math.sqrt(perPoint.reduce((a, b) => a + (b - mean) * (b - mean), 0) / (perPoint.length - 1))
        : NaN,
      relErr: 100 * (muB - P.MU_B) / P.MU_B,
      muBCodata: P.MU_B,
      n, ref: this.ref
    };
  };

  /** Nearest detected ring to a target radius, optionally within a tolerance. */
  function nearestRing(rings, target, tol) {
    let best = null, bd = Infinity;
    rings.forEach(p => {
      const d = Math.abs(p.r - target);
      if (d < bd) { bd = d; best = p; }
    });
    if (!best) return null;
    if (tol != null && bd > tol) return null;
    return best;
  }

  /**
   * All rings of one interference order, innermost first, tagged with their
   * 1..9 line index.  Used by the renderer to label the pattern so a student
   * can see which ring is "line 4" and which is "line 6".
   */
  function orderRings(m, B, fp, thetaMax) {
    const out = [];
    for (let i = 1; i <= 9; i++) {
      const r = P.splitDiameter(m, i, B, fp) / 2;
      out.push({ idx: i, r, delta: P.deltaOfLineIndex(i) });
    }
    return out.sort((a, b) => a.r - b.r);
  }

  root.ZeemanMeasure = { Session, orderRings };
})(typeof self !== 'undefined' ? self : this);
