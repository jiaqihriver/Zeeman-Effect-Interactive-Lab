/* Headless verification of the Zeeman simulator.
   Checks: no JS errors, physics values correct, measurement workflow recovers
   mu_B, and captures screenshots. */
const puppeteer = require('puppeteer-core');
const path = require('path');
const fs = require('fs');

const CHROME = '/Applications/Google Chrome.app/Contents/MacOS/Google Chrome';
const FILE = 'file://' + path.resolve(__dirname, 'index.html');
const OUT = path.resolve(__dirname, 'shots');

(async () => {
  fs.mkdirSync(OUT, { recursive: true });
  const browser = await puppeteer.launch({
    executablePath: CHROME, headless: 'new',
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
  R.errors = errors;
  R.dialogs = dialogs;
  R.booted = await page.evaluate(() => !!window.__ZEEMAN__);

  if (R.booted) {
    R.physics = await page.evaluate(() => {
      const { P } = window.__ZEEMAN__;
      const c0 = P.components(0), c1 = P.components(1.0);
      return {
        lambda0: +P.LAMBDA0_NM.toFixed(4), nu0: +P.NU0.toFixed(3),
        gUp: P.G_UP, gLo: P.G_LO, lorentz: +P.LORENTZ_PER_T.toFixed(5),
        nComp: c0.length,
        sumStrength: +c0.reduce((a, b) => a + b.strength, 0).toFixed(6),
        deltas: c0.map(x => x.delta),
        span1T_cm1: +(4 * P.LORENTZ_PER_T).toFixed(4),
        muB: P.MU_B
      };
    });

    // ---------- 2. canvas actually painted (not blank) ----------
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

    // ---------- 3. rings change with B (causal chain) ----------
    const ringHash = async (B) => page.evaluate((b) => {
      const { st, setB } = window.__ZEEMAN__;
      setB(b);
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

    // ---------- 4. geometry gating ----------
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

    // ---------- 5. two-step measurement workflow ----------
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

    // ---------- 6. teacher mode reveals result ----------
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
    // student mode must hide the answer
    await page.evaluate(() => document.querySelector('#segMode button[data-m="student"]').click());
    await new Promise(r => setTimeout(r, 200));
    R.studentModeHides = await page.evaluate(() =>
      !document.getElementById('resultBox').textContent.includes('CODATA'));

    // ---------- screenshots ----------
    const shot = async (name) => { await page.screenshot({ path: path.join(OUT, name), fullPage: true }); };
    await page.evaluate(() => {
      document.querySelector('#segMode button[data-m="teacher"]').click();
      window.__ZEEMAN__.setB(0.75);
    });
    await new Promise(r => setTimeout(r, 400));
    await shot('teacher-B0.75.png');

    // observation along B: pi lines vanish, only circular sigma remain
    await page.evaluate(() => {
      document.querySelector('#segView button[data-v="para"]').click();
      window.__ZEEMAN__.setB(0.9);
    });
    await new Promise(r => setTimeout(r, 400));
    await shot('teacher-along-B.png');
    await page.evaluate(() => {
      document.querySelector('#segView button[data-v="perp"]').click();
    });

    // B = 0 : nothing split
    await page.evaluate(() => {
      document.querySelector('#segMode button[data-m="student"]').click();
      window.__ZEEMAN__.setB(0);
    });
    await new Promise(r => setTimeout(r, 400));
    await shot('student-B0.png');
  }

  console.log(JSON.stringify(R, null, 1));
  await browser.close();
})();
