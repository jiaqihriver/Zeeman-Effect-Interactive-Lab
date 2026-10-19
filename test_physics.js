/* ==========================================================================
   test_physics.js — physics invariants for the Zeeman simulator.

   Runs on plain node with NO dependencies and NO browser:

       node test_physics.js          run the suite
       node test_physics.js --verbose  also print every passing check

   Why this file exists
   --------------------
   Every invariant below corresponds to a mistake this project has actually
   made, or to a convention that two modules previously disagreed about:

     - splitDiameter() used to recompute the wavenumber with the OPPOSITE sign
       to components(), mirroring the whole 1..9 ring labelling;
     - the field slider advertised a maximum the magnet model cannot reach;
     - analyse() converted every point with the currently selected pair.

   verify_headless.js covers rendering and the full measurement workflow, but
   it needs puppeteer and a specific Chrome install.  These checks must not
   depend on a browser being present, so the physics floor is enforced here.

   Note on loading: physics.js is a UMD module (module.exports under node),
   while measure.js reaches for the GLOBAL ZeemanPhysics.  Both are wired up
   below so the same source runs unchanged in the browser.
   ========================================================================== */
'use strict';

const path = require('path');
const g = globalThis;

g.self = g;
const P = require(path.join(__dirname, 'src', 'physics.js'));
g.ZeemanPhysics = P;
/* render.js needs a DOM; measure.js only needs these two pure functions, so
   stub them. Anything that genuinely needs pixels belongs in verify_headless. */
g.ZeemanRender = {
  radialProfile: () => ({ prof: new Float64Array(1), N: 1, rMax: 1 }),
  findRings: () => []
};
require(path.join(__dirname, 'src', 'measure.js'));
const M = g.ZeemanMeasure;

const VERBOSE = process.argv.includes('--verbose');
const TOL = 1e-9;

let passed = 0;
const failures = [];
let group = '';

function describe(name) { group = name; }
function check(name, condition, detail) {
  const label = group ? group + ' › ' + name : name;
  if (condition) {
    passed++;
    if (VERBOSE) console.log('  ok   ' + label + (detail ? '   ' + detail : ''));
  } else {
    failures.push(label + (detail ? '   [' + detail + ']' : ''));
    console.log('  FAIL ' + label + (detail ? '   ' + detail : ''));
  }
}
function near(name, got, want, tol, unit) {
  const d = Math.abs(got - want);
  check(name, d <= tol,
    'got ' + got + (unit ? ' ' + unit : '') + ', want ' + want +
    ' (Δ' + d.toExponential(1) + ')');
}

/* ------------------------------------------------------------- constants */
describe('atomic data (NIST ASD / CODATA 2018)');
near('lambda0 from level energies', P.LAMBDA0_NM, 546.227, 1e-2, 'nm');
near('NU0 = E_upper - E_lower', P.NU0, 18307.416, 1e-3, 'cm^-1');
near('Landé g(3S1)', P.G_UP, 2.0, TOL);
near('Landé g(3P2)', P.G_LO, 1.5, TOL);
near('Lorentz unit = μ_B B / hc', P.LORENTZ_PER_T, 0.4668645, 1e-6, 'cm^-1/T');
near('mu_B', P.MU_B, 9.2740100783e-24, 1e-33, 'J/T');

describe('Clebsch-Gordan strengths');
{
  const c = P.components(1.0);
  /* For fixed m_u the allowed m_l values form one column of the CG matrix,
     so the column must be normalised.  cg_exact.py derives the numbers and
     asserts the same invariant independently. */
  for (let mu = 1; mu >= -1; mu--) {
    const s = c.filter(x => x.mU === mu).reduce((a, b) => a + b.strength, 0);
    near('Σ|CG|² over fixed m_u = ' + mu, s, 1, 1e-12);
  }
  near('Σ all |CG|² = 2·J_u+1', c.reduce((a, b) => a + b.strength, 0), 3, 1e-12);
  check('exactly 9 allowed transitions', c.length === 9, 'got ' + c.length);
  check('every m_u slice is non-empty',
    [1, 0, -1].every(mu => c.some(x => x.mU === mu)));
}

/* ------------------------------------------------------ multiplet physics */
describe('multiplet structure');
{
  const c = P.components(1.0);
  const deltas = c.map(x => x.delta).sort((a, b) => a - b);
  check('9 distinct shifts',
    new Set(deltas.map(d => d.toFixed(6))).size === 9, JSON.stringify(deltas));
  check('shifts symmetric about 0',
    deltas.every((d, i) => Math.abs(d + deltas[8 - i]) < TOL));
  const sorted = c.slice().sort((a, b) => a.delta - b.delta);
  check('ν increases monotonically with δ',
    sorted.every((x, i) => i === 0 || x.nu_cm > sorted[i - 1].nu_cm));
  check('polarisation label matches Δm',
    c.every(x => (x.q === 0 && x.pol === 'pi') ||
                  (x.q > 0 && x.pol === 'sigma+') ||
                  (x.q < 0 && x.pol === 'sigma-')));

  const z0 = P.components(0).find(x => x.delta === 0);
  const z1 = P.components(1).find(x => x.delta === 0);
  near('δ = 0 line is field-independent', z1.nu_cm, z0.nu_cm, 1e-12, 'cm^-1');

  const at0 = P.components(0).map(x => x.lambda_nm);
  check('B = 0 is exactly degenerate (9 lines coincide)',
    new Set(at0.map(v => v.toFixed(12))).size === 1);

  const tiny = P.components(1e-9);
  const dev = Math.max(...tiny.map((x, i) => Math.abs(x.lambda_nm - at0[i])));
  check('splitting → 0 continuously as B → 0', dev < 1e-6, dev.toExponential(2) + ' nm');
}

/* ------------------------------------------------- ring order consistency */
describe('ring ordering: index ↔ shift ↔ radius  (regression)');
{
  /* This is the check that would have caught the mirrored-ring bug.
     components() derives ν from the level energies and is the authority on
     which wavelength belongs to which δ; splitDiameter() must agree with it
     to the last digit instead of recomputing ν on its own. */
  const m = 3657, B = 1.0;
  let worst = 0, worstAt = '';
  for (let i = 1; i <= 9; i++) {
    const d = P.deltaOfLineIndex(i);
    const fromIndex = P.splitDiameter(m, i, B, P.FP) / 2;
    const ref = P.components(B).find(x => Math.abs(x.delta - d) < 1e-9);
    if (!ref) { check('δ = ' + d + ' exists in components()', false); continue; }
    const fromComponents = P.ringRadius(m, ref.lambda_nm, P.FP);
    const dev = Math.abs(fromIndex - fromComponents);
    if (dev > worst) { worst = dev; worstAt = 'index ' + i + ' (δ = ' + d + ')'; }
  }
  check('splitDiameter agrees with components()', worst < TOL,
    worstAt + ' deviates by ' + worst.toExponential(2) + ' mm');

  /* Index 1 is the innermost ring. λ is longest when δ is most negative, and
     cos θ = mλ/2nd makes a longer λ a SMALLER angle, hence a smaller radius —
     so the reddest component must be innermost, not outermost. */
  const radii = [];
  for (let i = 1; i <= 9; i++) {
    radii.push({ i, d: P.deltaOfLineIndex(i), r: P.splitDiameter(m, i, B, P.FP) / 2 });
  }
  const byR = radii.slice().sort((a, b) => a.r - b.r);
  near('innermost ring (index 1) carries δ = −2', byR[0].d, -2, TOL);
  near('outermost ring (index 9) carries δ = +2', byR[8].d, 2, TOL);
  check('radius increases monotonically with index',
    radii.every((x, i) => i === 0 || x.r > radii[i - 1].r));

  /* The reddest line must have the longest wavelength and the smallest radius. */
  const c = P.components(B);
  const red = c.find(x => x.delta === -2), blue = c.find(x => x.delta === 2);
  check('δ = −2 has the longest wavelength', red.lambda_nm > blue.lambda_nm,
    red.lambda_nm.toFixed(4) + ' nm vs ' + blue.lambda_nm.toFixed(4) + ' nm');
  check('δ = −2 has the smallest ring radius',
    P.ringRadius(m, red.lambda_nm, P.FP) < P.ringRadius(m, blue.lambda_nm, P.FP));

  near('deltaOfLineIndex(5) is the unshifted line', P.deltaOfLineIndex(5), 0, TOL);
  P.LINE_INDEX_DELTA.forEach((d, i) => {
    near('lineIndexOfDelta round-trips δ = ' + d, P.lineIndexOfDelta(d), i + 1, 0);
  });
  check('out-of-range index yields undefined',
    P.deltaOfLineIndex(0) === undefined && P.deltaOfLineIndex(10) === undefined);
  check('unknown δ maps to null', P.lineIndexOfDelta(0.37) === null);
}

/* ------------------------------------------------------- measurement chain */
describe('μ_B extraction');
{
  const m = 3657, B = 1.0;
  P.SYMMETRIC_PAIRS.forEach(pr => {
    const Di = P.splitDiameter(m, pr.inner, B, P.FP);   // smaller ring
    const Do = P.splitDiameter(m, pr.outer, B, P.FP);   // larger ring
    check(pr.tag + ': inner ring is the smaller one', Di < Do,
      Di.toFixed(4) + ' mm vs ' + Do.toFixed(4) + ' mm');
    near(pr.tag + ': pair is symmetric about δ = 0',
      P.deltaOfLineIndex(pr.inner) + P.deltaOfLineIndex(pr.outer), 0, TOL);
    near(pr.tag + ': |δ| equals d0',
      Math.abs(P.deltaOfLineIndex(pr.inner)), pr.d0, TOL);
    const mu = P.muBFromDiameterPair(Do, Di, pr.d0, m, B, P.FP);
    near(pr.tag + ': recovers CODATA μ_B', mu, P.MU_B, 1e-6 * P.MU_B, 'J/T');
  });

  check('guarded against non-positive B, m, d0',
    Number.isNaN(P.muBFromDiameterPair(27, 26, 1, m, 0, P.FP)) &&
    Number.isNaN(P.muBFromDiameterPair(27, 26, 1, 0, B, P.FP)) &&
    Number.isNaN(P.muBFromDiameterPair(27, 26, 0, m, B, P.FP)));

  /* Splitting must grow linearly in B, so the linear fit the lab performs is
     legitimate. Check DSq at B and 2B doubles. */
  const pr = P.SYMMETRIC_PAIRS[1];
  const dsq = b => {
    const i = P.splitDiameter(m, pr.inner, b, P.FP);
    const o = P.splitDiameter(m, pr.outer, b, P.FP);
    return o * o - i * i;
  };
  near('ΔD² is linear in B (DSq(2B) = 2·DSq(B))', dsq(2 * B), 2 * dsq(B),
    1e-6 * Math.abs(2 * dsq(B)), 'mm²');
}

/* --------------------------------------------------- measurement session */
describe('measurement session');
{
  const mk = () => {
    const s = new M.Session();
    s.noiseMm = 0;                      // isolate the model from read noise
    s.geom = { view: 'perp', polariser: null, qwp: null };
    s.ref = { k: 3658, u: 1e-7, uMeasured: 1e-7, uRelError: 0, kExact: 3658,
              cosK: 1, cosK1: 1, residual: 0, Dk: 1, Dk1: 1, source: 'test' };
    return s;
  };
  const row = (B, d0) => ({ B, m: 3657, pairIdx: 0, d0, DSq: 1, muB: P.MU_B,
                            I: 0, Douter: 0, Dinner: 0, rOuter: 0, rInner: 0,
                            driftInner: 0, driftOuter: 0,
                            sepMeasured: 0, sepTheory: 0 });

  const s = mk();
  s.setPair(1);
  s.rows = [row(0.5, 1), row(1.0, 1)];
  check('analyse() accepts a consistent data set', !!s.analyse() && !s.analyse().error);

  /* Regression: rows measured with one pair must not be converted with
     another pair's prefactor — that scales μ_B by exactly the d0 ratio. */
  const mixed = mk();
  mixed.rows = [row(0.5, 0.5), row(1.0, 0.5)];   // d0 = 0.5 data ...
  check('analyse() rejects rows from a different pair',
    !!(mixed.analyse() && mixed.analyse().error),
    'd0 = 0.5 data with the d0 = 1 pair selected');

  check('setPair() drops stale rows',
    (() => { const t = mk(); t.rows = [row(1, 1)]; t.setPair(0); return t.rows.length === 0; })());

  check('analyse() needs a reference', new M.Session().analyse() === null);
  check('analyse() needs at least two points',
    (() => { const t = mk(); t.rows = [row(1, 1)]; return t.analyse() === null; })());
  check('session defaults to a seeded, reproducible PRNG',
    new M.Session().seed === 20261004 && new M.Session().noiseMm > 0);
}

/* ------------------------------------------------------------- geometry */
describe('observation geometry and polarisation');
{
  const c = P.components(1.0);
  const count = (view, pol) =>
    P.applyGeometry(c, view, pol).filter(x => x.weight >= P.VISIBLE_EPS).length;

  /* Signature is polariserWeight(polariser, isPi): the ANGLE comes first and
   isPi (the π-line flag) second. Viewed ⊥B, π is polarised ∥B and σ ⊥B, so
   Malus's law gives π → cos²φ and σ → sin²φ. A polariser aligned with B
   therefore passes the π lines and blocks σ, and vice versa at 90°. */
  near('no polariser ⇒ full transmission', P.polariserWeight(null, true), 1, TOL);
  near('polariser ∥B passes π fully', P.polariserWeight(0, true), 1, TOL);
  near('polariser ∥B kills σ', P.polariserWeight(0, false), 0, 1e-12);
  near('polariser ⊥B kills π', P.polariserWeight(90, true), 0, 1e-12);
  near('polariser ⊥B passes σ fully', P.polariserWeight(90, false), 1, TOL);
  near('π transmission varies continuously with angle',
    P.polariserWeight(45, true), 0.5, 1e-12);

  /* The bug normalisePolariser() exists to prevent: 0 is falsy, so a naive
     `opts.polariser || 'none'` silently shows the unfiltered multiplet. */
  check('polariser angle 0 is honoured, not treated as "none"',
    count('perp', 0) < count('perp', null) && P.normalisePolariser(0) === 0);
  check('"none" maps to null', P.normalisePolariser('none') === null);
  near('"parallel" maps to 0°', P.normalisePolariser('parallel'), 0, TOL);
  near('"perpendicular" maps to 90°', P.normalisePolariser('perpendicular'), 90, TOL);
  near('angle wraps into 0..360', P.normalisePolariser(370), 10, TOL);
  near('negative angle wraps', P.normalisePolariser(-90), 270, TOL);

  check('viewed along B, π lines do not radiate',
    P.applyGeometry(c, 'para').filter(x => x.q === 0)
      .every(x => x.weight === 0));
  near('linear polariser halves circular light along B',
    P.applyGeometry(c, 'para', 30).filter(x => x.q !== 0)[0].weight, 0.5, 1e-12);
  check('no component carries weight outside 0..1',
    P.applyGeometry(c, 'perp', 37).every(x => x.weight >= 0 && x.weight <= 1));

  /* QWP as a circular-polarisation separator. */
  near('QWP alone passes everything', P.qwpWeight(20, null, 1), 1, TOL);
  near('QWP + polariser at pol/2 passes σ+ and blocks σ−',
    P.qwpWeight(30, 60, 1), 1, 1e-9);
  near('QWP + polariser at pol/2 blocks σ+',
    P.qwpWeight(30, 60, -1), 0, 1e-9);
  check('qwpIsSeparator() recognises the separating setting',
    P.qwpIsSeparator(30, 60) && !P.qwpIsSeparator(22, 60));
  check('QWP is ignored when not inserted', P.qwpWeight(null, 60, 1) === null);
  near('QWP angle wraps into 0..180', P.normaliseQWP(190), 10, TOL);
}

/* ----------------------------------------------------------- instrument */
describe('interferometer and magnet');
{
  const fp = P.FP;
  near('finesse from R', fp.finesse, Math.PI * Math.sqrt(fp.R) / (1 - fp.R), 1e-9);
  near('Airy F from R', fp.F, 4 * fp.R / ((1 - fp.R) ** 2), 1e-9);
  near('FSR = 1/(2 n d)', fp.fsr_cm, 1 / (2 * fp.n * (fp.d_mm / 10)), 1e-12);

  /* Airy transmission must peak exactly at the resonance order m used by
     orderAt()/ringRadius(); the factor of 2 is easy to get wrong. */
  const lam = P.LAMBDA0_NM;
  let bestErr = 0;
  for (let m = 3654; m <= 3658; m++) {
    const theta = Math.acos(m / (2 * fp.n * (fp.d_mm / 10) / (lam * 1e-7)));
    const t = P.airy(Math.cos(theta), lam, fp);
    bestErr = Math.max(bestErr, Math.abs(t - 1));
  }
  check('Airy peaks at the resonance order', bestErr < 1e-9,
    'max |T − 1| = ' + bestErr.toExponential(2));

  near('orderAt inverts ringRadius',
    P.ringRadius(P.orderAt(0.01, lam, fp), lam, fp),
    fp.f_mm * Math.tan(0.01), 1e-9, 'mm');
  check('ringRadius rejects non-physical cos θ',
    Number.isNaN(P.ringRadius(1e9, lam, fp)));

  near('maxOrder is the on-axis order', P.maxOrder(fp),
    Math.floor(2 * fp.n * (fp.d_mm / 10) / (lam * 1e-7)), 1e-9);
  check('minVisibleOrder ≤ maxOrder',
    P.minVisibleOrder(4 * Math.PI / 180, fp) <= P.maxOrder(fp));

  /* The B slider must not advertise fields the magnet cannot produce:
     st.B = 1.30 tanh(I/5.20) saturates, so I_max caps the reachable field. */
  near('max reachable B at I_max', P.fieldFromCurrent(P.MAGNET.I_max), 1.0647, 1e-3, 'T');
  check('field is monotonic in current',
    P.fieldFromCurrent(1) < P.fieldFromCurrent(3) &&
    P.fieldFromCurrent(3) < P.fieldFromCurrent(6));
  check('negative current is clamped to zero field', P.fieldFromCurrent(-2) === 0);
  /* B = 1.05 T (the slider maximum in index.html) must be reachable at I ≤ 6 A. */
  const r = Math.min(0.999, 1.05 / P.MAGNET.B_sat);
  check('slider maximum B = 1.05 T is reachable within I_max',
    P.MAGNET.I0 * Math.atanh(r) <= P.MAGNET.I_max,
    'needs ' + (P.MAGNET.I0 * Math.atanh(r)).toFixed(2) + ' A of ' +
    P.MAGNET.I_max + ' A');
}

describe('linear fit');
{
  const xs = [0.5, 0.8, 1.0, 1.2], ys = xs.map(x => 3 * x + 1);
  const f = P.linfit(xs, ys);
  near('slope', f.slope, 3, 1e-9);
  near('intercept', f.intercept, 1, 1e-9);
  near('r²', f.r2, 1, 1e-12);
  check('rejects fewer than two points', P.linfit([1], [1]) === null);
  check('rejects degenerate input (all x equal)', P.linfit([1, 1, 1], [1, 2, 3]) === null);
  near('perfect fit has ~zero slope error', f.slopeErr, 0, 1e-9);
}

/* ------------------------------------------------------------------ done */
console.log('\n' + '-'.repeat(66));
if (failures.length) {
  console.log('FAILED  ' + failures.length + ' of ' + (passed + failures.length) +
              ' checks (' + passed + ' passed)');
  failures.forEach(f => console.log('   · ' + f));
  console.log('-'.repeat(66));
  process.exit(1);
}
console.log('PASSED  all ' + passed + ' checks');
console.log('-'.repeat(66));