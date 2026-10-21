/* Headless verification of the Zeeman simulator.
   Acts as a CI GATE: every check below is asserted, and the process exits
   with a non-zero status if any of them fails.  (The previous version only
   printed its findings, so a completely broken page still reported success.)

   Checks: no JS errors, physics values correct, the line-index labels sit on
   the rings they annotate, the measurement protocol is enforced, the
   measurement workflow recovers mu_B, and the dual-mode separation works.
   Screenshots are written to shots/ (gitignored). */
const puppeteer = require('puppeteer-core');
const path = require('path');
const fs = require('fs');

/* Locate a Chrome/Chromium binary.  The previous version hard-coded one
   macOS path, which made the script unusable on any other machine -- and
   therefore useless as a CI gate. */
const CHROME_CANDIDATES = [
  process.env.CHROME_PATH,
  '/Applications/Google Chrome.app/Contents/MacOS/Google Chrome',
  '/Applications/Chromium.app/Contents/MacOS/Chromium',
  '/usr/bin/google-chrome',
  '/usr/bin/google-chrome-stable',
  '/usr/bin/chromium',
  '/usr/bin/chromium-browser',
  'C:\\Program Files\\Google\\Chrome\\Application\\chrome.exe',
  'C:\\Program Files (x86)\\Google\\Chrome\\Application\\chrome.exe'
].filter(Boolean);

function findChrome() {
  for (const p of CHROME_CANDIDATES) {
    try { if (fs.existsSync(p)) return p; } catch (e) { /* ignore */ }
  }
  return null;
}

const FILE = 'file://' + path.resolve(__dirname, 'index.html');
const OUT = path.resolve(__dirname, 'shots');

/* ------------------------------------------------------------------ asserts */
const failures = [];
const checks = [];
function check(name, cond, detail) {
  const ok = !!cond;
  checks.push({ name, ok, detail: detail === undefined ? '' : String(detail) });
  if (!ok) failures.push(name + (detail !== undefined ? '  [' + detail + ']' : ''));
}
function near(a, b, tol) { return Math.abs(a - b) <= tol; }

(async () => {
  const chrome = findChrome();
  if (!chrome) {
    console.error('FATAL: no Chrome/Chromium binary found.');
    console.error('Set CHROME_PATH to point at one, e.g.');
    console.error('  CHROME_PATH=/usr/bin/chromium npm test');
    process.exit(2);
  }

  fs.mkdirSync(OUT, { recursive: true });
  const browser = await puppeteer.launch({
    executablePath: chrome,
    headless: true,                    // 'new' is deprecated in puppeteer >= 22
    args: ['--allow-file-access-from-files', '--no-sandbox', '--force-device-scale-factor=1']
  });
  const page = await browser.newPage();
  await page.setViewport({ width: 1680, height: 1150, deviceScaleFactor: 1.4 });

  const errors = [], logs = [];
  page.on('pageerror', e => errors.push('PAGEERROR: ' + e.message));
  page.on('console', m => { if (m.type() === 'error') errors.push('CONSOLE: ' + m.text()); else logs.push(m.text()); });
  page.on('requestfailed', r => errors.push('REQFAIL: ' + r.url()));
  // Headless Chrome SUSPENDS on alert()/confirm() until a handler acknowledges
  // it, so an unattended run would hang forever.  Collect and dismiss.
  const dialogs = [];
  page.on('dialog', async d => { dialogs.push(d.message()); try { await d.dismiss(); } catch (e) {} });

  await page.goto(FILE, { waitUntil: 'networkidle0' });
  await new Promise(r => setTimeout(r, 700));

  const R = {};

  // ---------- 1. static checks ----------
  R.booted = await page.evaluate(() => !!window.__ZEEMAN__);
  check('page boots (window.__ZEEMAN__ ready)', R.booted);
  check('no JavaScript errors', errors.length === 0, errors.join(' | '));
  check('no blocking dialogs', dialogs.length === 0, dialogs.join(' | '));

  if (R.booted) {
    // ---------- 2. physics values ----------
    R.physics = await page.evaluate(() => {
      const { P } = window.__ZEEMAN__;
      const c0 = P.components(0);
      const sorted = c0.slice().sort((a, b) => b.delta - a.delta);
      return {
        lambda0: +P.LAMBDA0_NM.toFixed(4), nu0: +P.NU0.toFixed(3),
        gUp: P.G_UP, gLo: P.G_LO, lorentz: +P.LORENTZ_PER_T.toFixed(5),
        nComp: c0.length,
        sumStrength: +c0.reduce((a, b) => a + b.strength, 0).toFixed(6),
        deltas: c0.map(x => x.delta),
        // Relative strengths per line, in the order delta = +2 ... -2.
        // These must follow the CG squares 1:3:6:3:4:3:6:3:1 -- a mere sum
        // check (which is 3 for ANY permutation) would not catch a table
        // whose entries had been mixed up.
        strengthRatios: sorted.map(x => +(x.strength * 10).toFixed(3)),
        span1T_cm1: +(4 * P.LORENTZ_PER_T).toFixed(4),
        muB: P.MU_B,
        // the magnet curve must be self-consistent: the field reachable at
        // I_max, and its inverse, have to agree with the slider bounds
        magnet: {
          B_max: +P.MAGNET.B_max.toFixed(6),
          I_max: +P.MAGNET.I_max.toFixed(4),
          B_at_Imax: +P.fieldFromCurrent(P.MAGNET.I_max).toFixed(6),
          I_at_Bmax: +P.currentFromField(P.MAGNET.B_max).toFixed(4)
        }
      };
    });

    check('9 Zeeman components', R.physics.nComp === 9, R.physics.nComp);
    check('lambda0 = 546.2267 nm', near(R.physics.lambda0, 546.2267, 1e-3), R.physics.lambda0);
    check('nu0 = 18307.416 cm^-1', near(R.physics.nu0, 18307.416, 1e-3), R.physics.nu0);
    check('g(upper) = 2.000', near(R.physics.gUp, 2, 1e-9), R.physics.gUp);
    check('g(lower) = 1.500', near(R.physics.gLo, 1.5, 1e-9), R.physics.gLo);
    check('1 Lorentz unit = 0.46686 cm^-1/T', near(R.physics.lorentz, 0.46686, 1e-5), R.physics.lorentz);
    check('sum of CG strengths = 3.000', near(R.physics.sumStrength, 3, 1e-6), R.physics.sumStrength);
    const wantDeltas = [2, 1.5, 1, 0.5, 0, -0.5, -1, -1.5, -2];
    check('component shifts are +2 .. -2', JSON.stringify(R.physics.deltas) === JSON.stringify(wantDeltas),
      JSON.stringify(R.physics.deltas));
    const wantRatios = [1, 3, 6, 3, 4, 3, 6, 3, 1];
    check('relative strengths are 1:3:6:3:4:3:6:3:1',
      JSON.stringify(R.physics.strengthRatios) === JSON.stringify(wantRatios),
      JSON.stringify(R.physics.strengthRatios));
    check('total span at 1 T = 1.8675 cm^-1', near(R.physics.span1T_cm1, 1.8675, 1e-4), R.physics.span1T_cm1);
    check('magnet curve: B(I_max) matches B_max',
      near(R.physics.magnet.B_at_Imax, R.physics.magnet.B_max, 1e-3),
      R.physics.magnet.B_at_Imax + ' vs ' + R.physics.magnet.B_max);
    check('magnet curve: inverse is consistent',
      near(R.physics.magnet.I_at_Bmax, R.physics.magnet.I_max, 0.011),
      R.physics.magnet.I_at_Bmax + ' vs ' + R.physics.magnet.I_max);

    // ---------- 2b. the public API surface stays self-consistent ----------
    // These functions are not used by the app itself, but they are part of
    // the documented developer interface, so a sign slip or a stale
    // convention in them would mislead anyone who used them.
    R.api = await page.evaluate(() => {
      const { P } = window.__ZEEMAN__;
      const m = 3657, B = 0.9, D = 4;
      const comps = P.components(B);
      const out = {};

      // the deprecated single-ring helpers must still invert to mu_B
      const slope = P.ringSlope(m, D, P.FP);
      out.slopeMuB = P.muBFromSlope(slope, m, D, P.FP);

      // ...and the ring they track must be the REDDEST one, i.e. the
      // INNERMOST ring of the group.  (The old comment claimed "outermost",
      // which is backwards: a longer wavelength gives a larger cos(theta)
      // and hence a smaller angle.)
      const r = Math.sqrt(P.ringRadiusSqVsB(m, B, D, P.FP));
      let best = null, bd = Infinity;
      comps.forEach(c => {
        const rr = P.ringRadius(m, c.lambda_nm, P.FP);
        const d = Math.abs(rr - r);
        if (d < bd) { bd = d; best = c; }
      });
      out.trackedDelta = best.delta;

      // delta -> index -> delta must round-trip for all nine lines
      out.roundTripOk = [2, 1.5, 1, 0.5, 0, -0.5, -1, -1.5, -2]
        .every(d => P.deltaOfLineIndex(P.lineIndexOfDelta(d)) === d);

      // meanSpacing must work on raw components AND respect visibility
      out.spacingRaw = P.meanSpacing(comps);
      out.spacingPiOnly = P.meanSpacing(P.applyGeometry(comps, 'perp', 90, null));

      // visibleComponents must agree with the explicit filter
      out.visibleCount = P.visibleComponents(comps, 'perp', null, null).length;
      out.filterCount = P.applyGeometry(comps, 'perp', null, null)
        .filter(x => x.weight >= P.VISIBLE_EPS).length;

      // suggestOrder must land on a usable reference order
      out.suggestOrder = P.suggestOrder(4 * Math.PI / 180);
      out.minVisibleOrder = P.minVisibleOrder(4 * Math.PI / 180);
      out.maxOrder = P.maxOrder();
      return out;
    });
    check('deprecated single-ring helpers still invert to mu_B',
      near(R.api.slopeMuB, R.physics.muB, R.physics.muB * 1e-3), R.api.slopeMuB);
    check('single-ring helper tracks the REDDEST (innermost) ring',
      R.api.trackedDelta === -2, R.api.trackedDelta);
    check('delta <-> line-index round-trips for all 9 lines', R.api.roundTripOk);
    check('meanSpacing works on raw components', R.api.spacingRaw > 0, R.api.spacingRaw);
    check('meanSpacing honours the visibility filter',
      R.api.spacingPiOnly > R.api.spacingRaw, R.api.spacingPiOnly + ' vs ' + R.api.spacingRaw);
    check('visibleComponents agrees with the explicit filter',
      R.api.visibleCount === R.api.filterCount, R.api.visibleCount + ' vs ' + R.api.filterCount);
    check('suggestOrder lands on a usable reference order',
      R.api.suggestOrder >= R.api.minVisibleOrder && R.api.suggestOrder < R.api.maxOrder,
      R.api.suggestOrder);

    // ---------- 3. the line-index labels sit on the right rings ----------
    // The nine rings of one order are labelled 1..9 from the innermost
    // outwards.  Each label must name the component that is ACTUALLY drawn
    // at that radius -- a sign slip here would put "+2" on the reddest ring
    // and mislead every student who reads the picture.
    R.labels = await page.evaluate(() => {
      const { P } = window.__ZEEMAN__;
      const m = 3657, B = 0.9;
      const comps = P.components(B);
      const rows = [];
      let mismatches = 0;
      for (let i = 1; i <= 9; i++) {
        const rLab = P.splitDiameter(m, i, B, P.FP) / 2;
        let best = null, bd = Infinity;
        comps.forEach(c => {
          const r = P.ringRadius(m, c.lambda_nm, P.FP);
          const d = Math.abs(r - rLab);
          if (d < bd) { bd = d; best = c; }
        });
        const label = P.deltaOfLineIndex(i);
        const ok = Math.abs(label - best.delta) < 1e-9;
        if (!ok) mismatches++;
        rows.push({ i, label, actual: best.delta, pol: best.pol, ok });
      }
      // and the inner ring must really be the reddest one
      const r1 = P.splitDiameter(m, 1, B, P.FP) / 2;
      const r9 = P.splitDiameter(m, 9, B, P.FP) / 2;
      return { mismatches, rows, innerIsSmaller: r1 < r9 };
    });
    check('all 9 ring labels match the drawn rings', R.labels.mismatches === 0,
      R.labels.rows.filter(r => !r.ok).map(r => `#${r.i} label ${r.label} vs ${r.actual}`).join(', '));
    check('line 1 is the innermost ring', R.labels.innerIsSmaller);

    // ---------- 4. canvas actually painted (not blank) ----------
    /* Read only a window from each canvas rather than the whole bitmap:
       getImageData on a full-resolution retina canvas costs tens of
       milliseconds per call, and a 200x200 sample is plenty to tell a painted
       canvas from an empty one. */
    R.canvasPainted = await page.evaluate(() => {
      const out = {};
      ['cvSpec', 'cvRing'].forEach(id => {
        const c = document.getElementById(id);
        const ctx = c.getContext('2d');
        const sw = Math.min(200, c.width), sh = Math.min(200, c.height);
        const d = ctx.getImageData(0, 0, sw, sh).data;
        let nonBg = 0, sum = 0, n = 0;
        for (let i = 0; i < d.length; i += 4 * 7) {
          const v = d[i] + d[i + 1] + d[i + 2];
          sum += v; n++;
          if (v > 40) nonBg++;
        }
        out[id] = { w: c.width, h: c.height, nonBgFrac: +(nonBg / n).toFixed(3), meanVal: +(sum / n / 3).toFixed(1) };
      });
      out.energySvg = document.getElementById('energyBox').innerHTML.length;
      out.apparatusSvg = document.getElementById('apparatus').innerHTML.length;
      out.chips = document.getElementById('chips').children.length;
      out.readoutCells = document.getElementById('readout').children.length;
      return out;
    });
    check('spectrum canvas is painted', R.canvasPainted.cvSpec.nonBgFrac > 0.01, R.canvasPainted.cvSpec.nonBgFrac);
    check('ring canvas is painted', R.canvasPainted.cvRing.nonBgFrac > 0.01, R.canvasPainted.cvRing.nonBgFrac);
    check('energy-level SVG is non-empty', R.canvasPainted.energySvg > 1000, R.canvasPainted.energySvg);
    check('apparatus SVG is non-empty', R.canvasPainted.apparatusSvg > 1000, R.canvasPainted.apparatusSvg);
    check('readout has cells', R.canvasPainted.readoutCells > 0, R.canvasPainted.readoutCells);

    // ---------- 5. rings change with B (causal chain) ----------
    // setB(..., now=true) forces a synchronous repaint: the canvas is read
    // back immediately below, so it must not wait for the next animation
    // frame.
    const ringHash = async (B) => page.evaluate((b) => {
      const { st, setB } = window.__ZEEMAN__;
      setB(b, false, true);
      const c = document.getElementById('cvRing');
      const ctx = c.getContext('2d');
      // sample two 160x160 windows on the horizontal diameter, where the
      // rings cross, so a change in their radii always changes the checksum
      let s = 0;
      [0.25, 0.6].forEach(fx => {
        const x = Math.floor(c.width * fx), y = Math.floor(c.height * 0.5) - 80;
        const d = ctx.getImageData(x, y, 160, 160).data;
        for (let i = 0; i < d.length; i += 4 * 5) s += d[i] * 3 + d[i + 1] * 5 + d[i + 2] * 7;
      });
      return s;
    }, B);
    R.ringHash = { b0: await ringHash(0), b05: await ringHash(0.5), b10: await ringHash(1.0) };
    R.ringHash.distinct = new Set(Object.values(R.ringHash)).size === 3;
    check('ring pattern changes with B', R.ringHash.distinct, JSON.stringify(R.ringHash));

    // ---------- 6. geometry gating ----------
    R.geometry = await page.evaluate(() => {
      const { P } = window.__ZEEMAN__;
      const c = P.components(1.0);
      const cnt = (v, p) => P.applyGeometry(c, v, p).filter(x => x.weight >= P.VISIBLE_EPS).length;
      return {
        perp_none: cnt('perp', 'none'), perp_par: cnt('perp', 'parallel'),
        perp_perp: cnt('perp', 'perpendicular'),
        para_none: cnt('para', 'none'), para_par: cnt('para', 'parallel')
      };
    });
    check('perp, no polariser: 9 visible', R.geometry.perp_none === 9, R.geometry.perp_none);
    check('perp, polariser ||B: 3 visible (pi)', R.geometry.perp_par === 3, R.geometry.perp_par);
    check('perp, polariser _|_B: 6 visible (sigma)', R.geometry.perp_perp === 6, R.geometry.perp_perp);
    check('along B: 6 visible (pi does not radiate)', R.geometry.para_none === 6, R.geometry.para_none);

    // ---------- 7. the measurement protocol is enforced ----------
    // The bench protocol requires: perpendicular to B, polariser removed, QWP
    // out.  Anything else breaks the 1..9 line numbering, so the measurement
    // buttons must be blocked and recorded points must be invalidated.
    R.protocol = await page.evaluate(async () => {
      const { st } = window.__ZEEMAN__;
      const ids = ['btnRef', 'btnMeasure', 'btnSweep'];
      const enabled = () => ids.map(id => !document.getElementById(id).disabled);
      const out = { initiallyEnabled: enabled() };
      document.querySelector('#polField .mini[data-a="90"]').click();
      await new Promise(r => setTimeout(r, 120));
      out.withPolariser = enabled();
      out.rowsAfterPolariser = st.session.rows.length;
      document.querySelector('#polField .mini[data-a="off"]').click();
      await new Promise(r => setTimeout(r, 120));
      out.afterRemoval = enabled();
      document.querySelector('#segView button[data-v="para"]').click();
      await new Promise(r => setTimeout(r, 120));
      out.alongB = enabled();
      document.querySelector('#segView button[data-v="perp"]').click();
      await new Promise(r => setTimeout(r, 120));
      out.backToPerp = enabled();
      return out;
    });
    check('measurement enabled in the protocol geometry',
      R.protocol.initiallyEnabled.every(Boolean), JSON.stringify(R.protocol.initiallyEnabled));
    check('measurement blocked with a polariser inserted',
      R.protocol.withPolariser.every(v => !v), JSON.stringify(R.protocol.withPolariser));
    check('measurement blocked when looking along B',
      R.protocol.alongB.every(v => !v), JSON.stringify(R.protocol.alongB));
    check('measurement re-enabled when the geometry is restored',
      R.protocol.backToPerp.every(Boolean), JSON.stringify(R.protocol.backToPerp));

    // ---------- 8. two-step measurement workflow ----------
    // Step 1 fixes the interference order from the two zero-field rings;
    // step 2 reads a symmetric split pair at each field.  Both symmetric
    // pairs are exercised because they have different sensitivities.
    R.measure = await page.evaluate(() => {
      const { st, M, P } = window.__ZEEMAN__;
      const out = { pairs: [] };
      for (let pi = 0; pi < P.SYMMETRIC_PAIRS.length; pi++) {
        const s = new M.Session();
        s.geom = { view: 'perp', polariser: null, qwp: null };
        s.setPair(pi);
        s.noiseMm = 0;                     // isolate the model from read noise
        const r0 = s.measureReference(st.orderK, { thetaMaxDeg: st.fovDeg, N: 8000 });
        if (!r0.ok) { out.pairs.push({ tag: s.pair().tag, fail: r0.reason }); continue; }
        const rows = [];
        for (const B of [0.50, 0.65, 0.80, 0.95, 1.05]) {
          const r = s.measurePair(B, 0, { thetaMaxDeg: st.fovDeg, N: 8000 });
          rows.push(r && r.ok
            ? { B, DSq: +r.row.DSq.toFixed(4), Do: +r.row.Douter.toFixed(4), Di: +r.row.Dinner.toFixed(4) }
            : { B, fail: r && r.reason });
        }
        const a = s.analyse();
        out.pairs.push({
          tag: s.pair().tag, d0: s.pair().d0,
          orderK: s.ref.k, orderM: s.splitOrder(),
          u: s.ref.u, orderResidual: s.ref.residual,
          okRows: rows.filter(x => !x.fail).length,
          muB: a ? a.muB : null, r2: a ? a.fit.r2 : null,
          slope: a ? a.fit.slope : null, relErrPct: a ? a.relErr : null
        });
      }
      // the sigma pair is the default selection; report it as the headline
      out.headline = out.pairs[1] || out.pairs[0];
      return out;
    });
    const sigma = R.measure.pairs[1] || {};
    check('sigma pair: 5/5 field points measured', sigma.okRows === 5, sigma.okRows);
    check('sigma pair: mu_B within 1% of CODATA',
      sigma.muB !== null && Math.abs(sigma.relErrPct) < 1, sigma.relErrPct);
    check('sigma pair: R^2 > 0.999', sigma.r2 !== null && sigma.r2 > 0.999, sigma.r2);

    // ---------- 9. degenerate inputs must not crash ----------
    // Two points at the SAME field make the design matrix singular; the old
    // code dereferenced a null fit and threw.
    R.degenerate = await page.evaluate(() => {
      const { st, M } = window.__ZEEMAN__;
      const s = new M.Session();
      s.noiseMm = 0;
      s.measureReference(st.orderK, { thetaMaxDeg: st.fovDeg, N: 8000 });
      s.measurePair(0.80, 0, { thetaMaxDeg: st.fovDeg, N: 8000 });
      s.measurePair(0.80, 0, { thetaMaxDeg: st.fovDeg, N: 8000 });
      let crash = null;
      try { const a = s.analyse(); crash = { ok: true, muB: a ? a.muB : null }; }
      catch (e) { crash = { ok: false, error: e.constructor.name + ': ' + e.message }; }
      return crash;
    });
    check('duplicate field points do not crash the fit', R.degenerate.ok === true,
      R.degenerate.error || '');

    // ---------- 10. teacher mode reveals result ----------
    // Step 1 must precede the sweep, because the sweep needs the order that
    // step 1 determines.
    await page.evaluate(() => {
      document.querySelector('#segMode button[data-m="teacher"]').click();
      document.getElementById('btnRef').click();
    });
    await new Promise(r => setTimeout(r, 400));
    await page.evaluate(() => document.getElementById('btnSweep').click());
    await new Promise(r => setTimeout(r, 600));
    R.teacherMode = await page.evaluate(() => {
      const txt = document.getElementById('resultBox').textContent.replace(/\s+/g, ' ');
      return {
        bodyClass: document.body.classList.contains('teacher'),
        resultLen: document.getElementById('resultBox').innerHTML.length,
        showsSlope: txt.includes('斜率'),
        // the value is rendered as "9.1615 × 10⁻²⁴" with superscript glyphs
        showsMuB: /μ_B/.test(txt) && /9\.\d{4}/.test(txt),
        showsCodata: txt.includes('CODATA'),
        showsRelErr: txt.includes('相对误差'),
        snippet: txt.slice(0, 220)
      };
    });
    check('teacher mode reveals the slope', R.teacherMode.showsSlope);
    check('teacher mode reveals mu_B', R.teacherMode.showsMuB, R.teacherMode.snippet);
    check('teacher mode reveals the CODATA value', R.teacherMode.showsCodata);
    check('teacher mode reveals the relative error', R.teacherMode.showsRelErr);
    // student mode must hide the answer
    await page.evaluate(() => document.querySelector('#segMode button[data-m="student"]').click());
    await new Promise(r => setTimeout(r, 200));
    R.studentModeHides = await page.evaluate(() =>
      !document.getElementById('resultBox').textContent.includes('CODATA'));
    check('student mode hides the answer', R.studentModeHides);

    // ---------- screenshots ----------
    const shot = async (name) => { await page.screenshot({ path: path.join(OUT, name), fullPage: true }); };
    await page.evaluate(() => {
      document.querySelector('#segMode button[data-m="teacher"]').click();
      window.__ZEEMAN__.setB(0.75, false, true);
    });
    await new Promise(r => setTimeout(r, 400));
    await shot('teacher-B0.75.png');

    // observation along B: pi lines vanish, only circular sigma remain
    await page.evaluate(() => {
      document.querySelector('#segView button[data-v="para"]').click();
      window.__ZEEMAN__.setB(0.9, false, true);
    });
    await new Promise(r => setTimeout(r, 400));
    await shot('teacher-along-B.png');
    await page.evaluate(() => {
      document.querySelector('#segView button[data-v="perp"]').click();
    });

    // B = 0 : nothing split
    await page.evaluate(() => {
      document.querySelector('#segMode button[data-m="student"]').click();
      window.__ZEEMAN__.setB(0, false, true);
    });
    await new Promise(r => setTimeout(r, 400));
    await shot('student-B0.png');
  }

  // ------------------------------------------------------------------ report
  R.checks = checks;
  R.failures = failures;
  R.errors = errors;
  R.dialogs = dialogs;
  R.chrome = chrome;

  await browser.close();

  const passed = checks.filter(c => c.ok).length;
  console.log('');
  console.log('='.repeat(66));
  console.log('  Zeeman simulator -- headless regression');
  console.log('='.repeat(66));
  for (const c of checks) {
    console.log(`  ${c.ok ? 'PASS' : 'FAIL'}  ${c.name}${c.ok || !c.detail ? '' : '   -> ' + c.detail}`);
  }
  console.log('-'.repeat(66));
  console.log(`  ${passed}/${checks.length} checks passed` +
    (failures.length ? `   ${failures.length} FAILED` : '   all green'));
  console.log('='.repeat(66));

  if (failures.length) {
    console.error('\nFAILED:\n' + failures.map(f => '  - ' + f).join('\n'));
    process.exit(1);
  }
  if (errors.length) {
    console.error('\nRuntime errors:\n' + errors.map(e => '  - ' + e).join('\n'));
    process.exit(1);
  }
  process.exit(0);
})().catch(e => {
  console.error('FATAL: ' + (e && e.stack ? e.stack : e));
  process.exit(2);
});
