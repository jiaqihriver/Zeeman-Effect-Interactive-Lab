/* ==========================================================================
   render.js — canvas / SVG drawing for the Zeeman simulator
   Depends on physics.js (global ZeemanPhysics).
   ========================================================================== */
(function (root) {
  'use strict';
  const P = root.ZeemanPhysics;

  /* --------------------------------------------------------------- palette */
  const C = {
    bg: '#0b1020', panel: '#121a33', grid: 'rgba(148,180,255,0.10)',
    axis: '#8ea3cc', text: '#dce6ff', dim: '#8ea3cc',
    up: '#ff9f43', lo: '#4fc3f7', bField: '#a78bfa',
    pi: '#4ade80', sp: '#60a5fa', sm: '#f472b6',
    ok: '#4ade80', warn: '#fbbf24', bad: '#f87171'
  };
  const polColor = { 'pi': C.pi, 'sigma+': C.sp, 'sigma-': C.sm };

  function setupCanvas(canvas) {
    const dpr = window.devicePixelRatio || 1;
    // Measure the element's LAYOUT box.  The CSS gives every canvas an explicit
    // size, so this is stable and independent of the width/height attributes
    // (which only define the backing-store resolution).
    let rect = canvas.getBoundingClientRect();
    let w = rect.width, h = rect.height;
    if (!(w > 0 && h > 0)) {
      // first paint before layout: fall back to the parent box
      const p = canvas.parentElement;
      if (p) {
        const pr = p.getBoundingClientRect();
        w = pr.width; h = pr.height;
      }
    }
    w = Math.max(1, Math.round(w || 300));
    h = Math.max(1, Math.round(h || 300));
    const bw = Math.round(w * dpr), bh = Math.round(h * dpr);
    if (canvas.width !== bw || canvas.height !== bh) {
      canvas.width = bw; canvas.height = bh;
    }
    const ctx = canvas.getContext('2d');
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    ctx.clearRect(0, 0, w, h);
    return { ctx, w, h, dpr, bw, bh };
  }

  /* =======================================================================
     1. Energy level diagram (SVG) — upper 3S1 and lower 3P2 with Zeeman
        sublevels, transitions drawn as arrows coloured by polarisation.
     ======================================================================= */
  function energySVG(B, opts) {
    opts = opts || {};
    const comps = P.components(B);
    const vis = P.applyGeometry(comps, opts.view || 'perp', P.normalisePolariser(opts.polariser), P.normaliseQWP(opts.qwp));
    const Ju = P.LEVELS[P.LINE.upper].J, Jl = P.LEVELS[P.LINE.lower].J;
    const gU = P.G_UP, gL = P.G_LO;
    const active = Math.abs(B) > 1e-6;

    /* ------------------------------------------------------------------ plan
       Classic energy-level diagram, read left to right:

         |  ΔE table  |   sub-level lines (long)   |  m_J  |
         |  (values) |                            |       |

       Transitions are drawn as VERTICAL arrows in a dedicated column between
       the two manifolds.  Because every arrow spans only its own pair of
       sub-levels, no arrow can overshoot another and the picture stays
       readable at any field.  A summary table on the right lists the nine
       allowed transitions, so nothing has to be crammed onto the arrows. */
    const W = 1000;
    const AX = 92;                              // energy axis
    const TBL_X = 112, TBL_W = 208;              // Delta-E table column
    const LX0 = TBL_X + TBL_W + 14;              // 334: level lines start
    const LX1 = W - 172;                         // 828: level lines end
    const M_X = LX1 + 10;                        // m_J labels
    const yTop = 132, yBot = 340;                // zero-field level positions
    const SPLIT = 52;                            // px, full-scale deflection

    const shiftPx = (deltaLorentz) => {
      const dnu = deltaLorentz * P.LORENTZ_PER_T * B;             // cm^-1
      const maxShift = 2 * gU * P.LORENTZ_PER_T * Math.max(B, 1e-9);
      return (dnu / maxShift) * SPLIT;
    };
    // half-spread of each manifold, used to keep the bands and labels apart
    const sprU = (gU * Ju) * (SPLIT / (2 * gU));   // = Ju*SPLIT/2
    const sprL = (gL * Jl) * (SPLIT / (2 * gU));

    /* Vertical budget is computed up front so the canvas is exactly tall
       enough for the summary table — no row is ever clipped.               */
    const ROW_H = 14;
    const legY = yBot + sprL + 52;               // legend baseline
    const tbY = legY + 26;                       // table header baseline
    const tbX = TBL_X - 8;                      // table left edge
    const nRows = vis.length || 1;
    const H = Math.ceil(tbY + 15 + nRows * ROW_H + 26);

    let s = `<svg viewBox="0 0 ${W} ${H}" xmlns="http://www.w3.org/2000/svg" class="svg-fig">`;
    s += `<defs>
      <marker id="ah" viewBox="0 0 10 10" refX="8.5" refY="5" markerWidth="5.5" markerHeight="5.5" orient="auto">
        <path d="M0,0 L10,5 L0,10 z" fill="context-stroke"/>
      </marker></defs>`;
    s += `<rect x="0" y="0" width="${W}" height="${H}" fill="${C.bg}" rx="10"/>`;

    /* ---- energy axis ---- */
    const axT = yTop - sprU - 26, axB = yBot + sprL + 30;
    s += `<line x1="${AX}" y1="${axT}" x2="${AX}" y2="${axB}" stroke="${C.axis}" stroke-width="1.5"/>`;
    s += `<polygon points="${AX - 5},${axT} ${AX + 5},${axT} ${AX},${axT - 10}" fill="${C.axis}"/>`;
    s += `<text transform="rotate(-90 ${AX - 13} ${(axT + axB) / 2})" x="${AX - 13}" y="${(axT + axB) / 2}"
           fill="${C.dim}" font-size="14" text-anchor="middle">能量 E</text>`;

    /* ---- Delta-E column -------------------------------------------------
       One single line per sub-level.  Stacking the value and the Lorentz-unit
       value on two lines needs ~27 px of leading, but adjacent sub-levels are
       only ~13 px apart at full field, so the two lines would inevitably
       overlap.  A single line keeps the entry inside its own 13 px slot.     */
    const dERow = (yc, m, g) => {
      if (!active) return '';
      const dE = g * m * P.MU_B * B;
      return `<text x="${TBL_X + TBL_W}" y="${(yc + 4).toFixed(1)}" fill="${C.dim}"
              font-size="11" text-anchor="end">ΔE = ${fmtSci(dE)} J　(${fmtSigned(g * m, 2)} μ_B B)</text>`;
    };

    /* ---- upper manifold 6s7s 3S1 (J = 1) ---- */
    s += `<text x="${LX0}" y="${yTop - sprU - 26}" fill="${C.up}" font-size="15" font-weight="600"
           >6s7s ³S₁</text>
          <text x="${LX0 + 96}" y="${yTop - sprU - 26}" fill="${C.dim}" font-size="12"
           >J = 1　L = 0　S = 1　g = ${gU.toFixed(3)}</text>`;
    s += `<line x1="${LX0}" y1="${yTop}" x2="${LX1}" y2="${yTop}"
           stroke="${C.up}" stroke-width="1" opacity="0.28" stroke-dasharray="4 4"/>`;

    const uPos = {};
    for (let m = Ju; m >= -Ju; m--) {
      const y = yTop - shiftPx(gU * m);
      uPos[m] = y;
      s += `<line x1="${LX0}" y1="${y}" x2="${LX1}" y2="${y}" stroke="${C.up}"
             stroke-width="${active ? 2.4 : 2}" opacity="${active ? 1 : 0.5}"/>`;
      if (active) s += `<text x="${M_X}" y="${y + 4.5}" fill="${C.up}" font-size="13">m = ${m > 0 ? '+' + m : m}</text>`;
      s += dERow(y, m, gU);
    }
    /* B = 0: the sub-levels are coincident, so their m values go in ONE row
       above the manifold instead of stacking on the same line. */
    if (!active) {
      const list = Array.from({ length: 2 * Ju + 1 }, (_, k) => k - Ju)
        .map(v => v > 0 ? '+' + v : String(v)).join('  ·  ');
      s += `<text x="${M_X}" y="${yTop - sprU - 10}" fill="${C.up}" font-size="12.5">m = ${list}</text>`;
    }

    /* ---- lower manifold 6s6p 3P2 (J = 2) ---- */
    s += `<text x="${LX0}" y="${yBot + sprL + 30}" fill="${C.lo}" font-size="15" font-weight="600"
           >6s6p ³P₂</text>
          <text x="${LX0 + 96}" y="${yBot + sprL + 30}" fill="${C.dim}" font-size="12"
           >J = 2　L = 1　S = 1　g = ${gL.toFixed(3)}</text>`;
    s += `<line x1="${LX0}" y1="${yBot}" x2="${LX1}" y2="${yBot}"
           stroke="${C.lo}" stroke-width="1" opacity="0.28" stroke-dasharray="4 4"/>`;

    const lPos = {};
    for (let m = Jl; m >= -Jl; m--) {
      const y = yBot - shiftPx(gL * m);
      lPos[m] = y;
      s += `<line x1="${LX0}" y1="${y}" x2="${LX1}" y2="${y}" stroke="${C.lo}"
             stroke-width="${active ? 2.4 : 2}" opacity="${active ? 1 : 0.5}"/>`;
      if (active) s += `<text x="${M_X}" y="${y + 4.5}" fill="${C.lo}" font-size="13">m = ${m > 0 ? '+' + m : m}</text>`;
      s += dERow(y, m, gL);
    }
    if (!active) {
      const list = Array.from({ length: 2 * Jl + 1 }, (_, k) => k - Jl)
        .map(v => v > 0 ? '+' + v : String(v)).join('  ·  ');
      s += `<text x="${M_X}" y="${yBot + sprL + 14}" fill="${C.lo}" font-size="12.5">m = ${list}</text>`;
    }

    /* ---- transitions: vertical arrows, one lane each --------------------
       Lanes are spread across the middle of the level-line band and each
       arrow connects exactly one pair of sub-levels, so arrows never cross
       and never run past the lines.                                      */
    const bandT = yTop + sprU * 0.25, bandB = yBot - sprL * 0.25;
    const lanes = vis.length || 1;
    const laneW = Math.min(58, (LX1 - LX0 - 40) / lanes);
    const x0 = LX0 + (LX1 - LX0 - laneW * lanes) / 2;   // centred block
    vis.forEach((c, i) => {
      if (c.weight < P.VISIBLE_EPS) return;
      const y1 = uPos[c.mU], y2 = lPos[c.mL];
      if (y1 === undefined || y2 === undefined) return;
      const x = x0 + laneW * (i + 0.5);
      const col = polColor[c.pol];
      const isPi = c.q === 0;
      // stop clear of the arrow head so it does not touch the target line
      const yA = y1 + (y1 < y2 ? 4 : -4);
      const yB = y2 + (y1 < y2 ? -5 : 5);
      s += `<line x1="${x.toFixed(1)}" y1="${yA.toFixed(1)}" x2="${x.toFixed(1)}" y2="${yB.toFixed(1)}"
             stroke="${col}" stroke-width="${isPi ? 2.2 : 1.6}" opacity="${(0.5 + 0.5 * c.weight).toFixed(2)}"
             marker-end="url(#ah)"/>`;
    });

    /* ---- legend + transition summary table ------------------------------
       Vertical budget below the lower manifold, top to bottom:
         legend row -> table header -> 9 rows -> footer
       All of it is anchored to the lower manifold, so growing SPLIT pushes the
       block down instead of letting its parts collide.                     */
    s += `<rect x="${TBL_X - 8}" y="${legY - 15}" width="${TBL_W + 16}" height="25" rx="6"
           fill="#0d1428" opacity="0.9" stroke="rgba(140,170,240,.14)"/>`;
    const leg = [
      { t: 'σ⁺', c: C.sp, d: 'Δm = +1' },
      { t: 'π', c: C.pi, d: 'Δm = 0' },
      { t: 'σ⁻', c: C.sm, d: 'Δm = −1' }
    ];
    leg.forEach((it, i) => {
      const lx = TBL_X + i * 70;
      s += `<circle cx="${lx + 5}" cy="${legY - 2.5}" r="4.5" fill="${it.c}"/>`;
      s += `<text x="${lx + 14}" y="${legY + 1.5}" fill="${C.dim}" font-size="11.5">${it.t} ${it.d}</text>`;
    });

    // relative strengths as a small bar strip, to the right of the legend
    const sxs = LX0 + 24, sbY = legY;
    s += `<text x="${sxs}" y="${sbY + 1.5}" fill="${C.dim}" font-size="11.5">相对强度</text>`;
    let bx = sxs + 58;
    comps.forEach(c => {
      const hgt = 12 * (c.strength / 0.6);
      s += `<rect x="${bx.toFixed(1)}" y="${(sbY + 3 - hgt).toFixed(1)}" width="8" height="${hgt.toFixed(1)}"
             fill="${polColor[c.pol]}" opacity="${c.weight >= P.VISIBLE_EPS ? 0.92 : 0.2}"/>`;
      bx += 11;
    });

    /* ---- transition summary table ---------------------------------------
       Full width, below the diagram.  The canvas is sized from this budget so
       that all nine rows are always drawn — never silently clipped.        */
    const sorted = vis.slice().sort((a, b) => a.delta - b.delta || a.mL - b.mL);
    const nVis = sorted.filter(c => c.weight >= P.VISIBLE_EPS).length;
    const colX = [0, 46, 146, 216, 262, 340];
    s += `<text x="${tbX + colX[0]}" y="${tbY}" fill="${C.dim}" font-size="10.5">Δ / L</text>`;
    s += `<text x="${tbX + colX[1]}" y="${tbY}" fill="${C.dim}" font-size="10.5">λ / nm</text>`;
    s += `<text x="${tbX + colX[2]}" y="${tbY}" fill="${C.dim}" font-size="10.5">m_J 上 → 下</text>`;
    s += `<text x="${tbX + colX[3]}" y="${tbY}" fill="${C.dim}" font-size="10.5">Δm</text>`;
    s += `<text x="${tbX + colX[4]}" y="${tbY}" fill="${C.dim}" font-size="10.5">偏振</text>`;
    s += `<text x="${tbX + colX[5]}" y="${tbY}" fill="${C.dim}" font-size="10.5">相对强度</text>`;
    sorted.forEach((c, i) => {
      const y = tbY + 15 + i * ROW_H;
      const on = c.weight >= P.VISIBLE_EPS;
      const col = on ? polColor[c.pol] : C.dim;
      const op = on ? 1 : 0.4;
      const sg = c.q === 0 ? 'π' : c.q > 0 ? 'σ⁺' : 'σ⁻';
      s += `<text x="${tbX + colX[0]}" y="${y.toFixed(1)}" fill="${C.dim}" font-size="10.5" opacity="${on ? 0.85 : 0.4}"`
         + `>${c.delta > 0 ? '+' : ''}${c.delta.toFixed(1)}</text>`;
      s += `<text x="${tbX + colX[1]}" y="${y.toFixed(1)}" fill="${col}" font-size="10.5" opacity="${op}"`
         + `>${c.lambda_nm.toFixed(3)}</text>`;
      s += `<text x="${tbX + colX[2]}" y="${y.toFixed(1)}" fill="${col}" font-size="10.5" opacity="${op}"`
         + `>${c.mU > 0 ? '+' + c.mU : c.mU} → ${c.mL > 0 ? '+' + c.mL : c.mL}</text>`;
      s += `<text x="${tbX + colX[3]}" y="${y.toFixed(1)}" fill="${col}" font-size="10.5" opacity="${op}"`
         + `>${c.q > 0 ? '+' + c.q : c.q}</text>`;
      s += `<text x="${tbX + colX[4]}" y="${y.toFixed(1)}" fill="${col}" font-size="10.5" opacity="${op}"`
         + `>${sg}${on ? '' : ' ✕'}</text>`;
      // strength as a small inline bar plus the percentage
      const bw = 34 * (c.strength / 0.6);
      s += `<rect x="${tbX + colX[5]}" y="${(y - 8).toFixed(1)}" width="34" height="7" rx="2"
             fill="#0d1428" opacity="0.8"/>`;
      s += `<rect x="${tbX + colX[5]}" y="${(y - 8).toFixed(1)}" width="${bw.toFixed(1)}" height="7" rx="2"
             fill="${polColor[c.pol]}" opacity="${on ? 0.9 : 0.25}"/>`;
      s += `<text x="${tbX + colX[5] + 40}" y="${y.toFixed(1)}" fill="${col}" font-size="10.5" opacity="${op}"`
         + `>${(c.strength * 100).toFixed(0)}%</text>`;
    });

    /* ---- footer --------------------------------------------------------- */
    const fY = H - 10;
    s += `<text x="${TBL_X - 8}" y="${fY}" fill="${C.text}" font-size="12"
           >λ₀ = ${P.LINE.lambda0_nm_air.toFixed(2)} nm　ν₀ = ${P.NU0.toFixed(2)} cm⁻¹`
         + `${active ? `　总跨度 Δν = ${(4 * P.LORENTZ_PER_T * B).toFixed(4)} cm⁻¹` : ''}</text>`;
    if (!active) {
      // Left column only (B = 0 leaves the Delta-E column empty), kept short
      // enough not to run into the manifold titles further right.
      s += `<text x="${TBL_X - 8}" y="${yTop - sprU - 24}" fill="${C.warn}" font-size="11.5"
             >B = 0：无外场，9 条谱线重合为 1 条</text>`;
    }
    s += `<text x="${LX1}" y="${fY}" fill="${C.dim}" font-size="11" text-anchor="end"
           >ΔE = g_J m_J μ_B B（纵向刻度放大约 10⁵ 倍）</text>`;
    s += `</svg>`;
    return s;
  }

  /* =======================================================================
     2. Spectrum canvas — the 9 components as intensity peaks, with an
        optional unresolved/blurred view controlled by instrument resolution.
     ======================================================================= */
  function spectrum(canvas, B, opts) {
    opts = opts || {};
    const { ctx, w, h } = setupCanvas(canvas);
    const comps = P.applyGeometry(P.components(B), opts.view || 'perp', P.normalisePolariser(opts.polariser), P.normaliseQWP(opts.qwp));

    ctx.fillStyle = C.bg; ctx.fillRect(0, 0, w, h);
    const ML = 66, MR = 18, MT = 46, MB = 60;
    const pw = w - ML - MR, ph = h - MT - MB;

    // span: +-5 Lorentz units around nu0 (or +-6 for headroom)
    const spanLorentz = opts.spanLorentz || 2.6;
    const half = spanLorentz * P.LORENTZ_PER_T * Math.max(B, 0.12);
    const nuMin = P.NU0 - half, nuMax = P.NU0 + half;

    const x = nu => ML + (nu - nuMin) / (nuMax - nuMin) * pw;
    // intensity scale: normalise to the strongest visible component
    let smax = 0;
    comps.forEach(c => { if (c.weight >= P.VISIBLE_EPS) smax = Math.max(smax, c.strength); });
    smax = smax || 1;
    const y = v => MT + ph - (v / smax) * ph * 0.88;

    // grid
    ctx.strokeStyle = C.grid; ctx.lineWidth = 1;
    ctx.font = '12px ui-monospace, Menlo, monospace';
    ctx.fillStyle = C.dim;
    const step = half > 0.6 ? 0.2 : half > 0.25 ? 0.1 : 0.05;
    for (let d = -Math.floor(spanLorentz / step) * step; d <= spanLorentz + 1e-9; d += step) {
      const nu = P.NU0 + d * P.LORENTZ_PER_T * Math.max(B, 0.12);
      const px = x(nu);
      if (px < ML - 1 || px > ML + pw + 1) continue;
      const major = Math.abs(d) < 1e-9;
      ctx.beginPath(); ctx.moveTo(px, MT); ctx.lineTo(px, MT + ph);
      ctx.strokeStyle = major ? 'rgba(148,180,255,0.32)' : C.grid;
      ctx.stroke();
      if (major || Math.abs(d / step) % 2 < 1e-6) {
        ctx.fillText((d >= 0 ? '+' : '') + d.toFixed(2), px - 13, h - 44);
      }
    }
    // zero line
    ctx.strokeStyle = C.axis; ctx.lineWidth = 1.4;
    ctx.beginPath(); ctx.moveTo(ML, y(0)); ctx.lineTo(ML + pw, y(0)); ctx.stroke();

    // axis labels
    ctx.fillStyle = C.text; ctx.font = '13px system-ui, sans-serif';
    ctx.fillText('相对强度', ML, 15);
    ctx.textAlign = 'right';
    ctx.fillText('波数偏移 ν − ν₀  (单位: 洛伦兹单位 L = μ_B B/hc)', ML + pw, h - 8);
    ctx.textAlign = 'left';

    // resolution blur in Lorentz units
    const resL = opts.resolutionLorentz || 0.045;

    // Draw each visible component as a broadened peak.
    comps.forEach(c => {
      // Threshold, not `weight <= 0`: cos^2(90 deg) is 3.7e-33 in floating point,
      // not exactly zero, so a `>0` test still draws a zero-height sliver.
      if (c.weight < P.VISIBLE_EPS) return;
      const amp = (c.strength / smax) * c.weight;
      const col = polColor[c.pol];
      const x0 = x(c.nu_cm);
      // integrate Gaussian over the pixel span
      ctx.beginPath();
      const halfWidthPx = Math.max(1.2, Math.abs(x(c.nu_cm + resL * P.LORENTZ_PER_T * Math.max(B, 0.12)) - x0));
      ctx.moveTo(x0 - halfWidthPx * 2.2, y(0));
      for (let t = -2.2; t <= 2.2; t += 0.05) {
        const g = Math.exp(-0.5 * Math.pow(t / 0.55, 2));
        ctx.lineTo(x0 + halfWidthPx * t, y(amp * g));
      }
      ctx.lineTo(x0 + halfWidthPx * 2.2, y(0));
      ctx.closePath();
      const grad = ctx.createLinearGradient(0, y(0), 0, y(amp));
      grad.addColorStop(0, col + '00'); grad.addColorStop(1, col + 'ff');
      ctx.fillStyle = grad; ctx.fill();
      ctx.strokeStyle = col; ctx.lineWidth = 1.6; ctx.stroke();

      // peak marker + label
      ctx.beginPath(); ctx.moveTo(x0, y(amp)); ctx.lineTo(x0, y(0));
      ctx.strokeStyle = col; ctx.lineWidth = 1; ctx.setLineDash([2, 3]); ctx.stroke();
      ctx.setLineDash([]);
      if (opts.labels !== false) {
        ctx.fillStyle = col; ctx.font = '12.5px system-ui, sans-serif';
        ctx.textAlign = 'center';
        const tag = c.pol === 'pi' ? 'π' : c.pol === 'sigma+' ? 'σ⁺' : 'σ⁻';
        ctx.fillText(tag, x0, y(amp) - 7);
        ctx.textAlign = 'left';
      }
    });

    /* ---- sigma-group annotation (only meaningful viewed along B) ---------
       Along B thepi lines do not radiate, so the multiplet reduces to two
       circularly-polarised groups of three.  Bracket each group and label it
       with the upper sub-levels that feed it, which is the "where does this
       group come from" answer a lab report needs. */
    if (opts.groupBrackets && (opts.view || 'perp') === 'para' && B > 1e-6) {
      const grp = P.sigmaGroups(comps);
      const brY = MT - 4;                // bracket rail just under the axis title
      const vis2 = grp['sigma+'].vis + grp['sigma-'].vis;
      if (vis2 > 0) {
        [['sigma-', grp['sigma-'], C.sm], ['sigma+', grp['sigma+'], C.sp]]
          .forEach(([name, g, colr]) => {
            if (!g.n) return;
            const live = g.lines.filter(c => c.weight >= P.VISIBLE_EPS);
            if (!live.length) return;
            const xs = live.map(c => x(c.nu_cm));
            const x0b = Math.min.apply(null, xs), x1b = Math.max.apply(null, xs);
            const dead = g.vis < g.n;
            ctx.save();
            ctx.globalAlpha = dead ? 0.45 : 1;
            // rail with end ticks
            ctx.strokeStyle = colr; ctx.lineWidth = 1.4;
            ctx.beginPath();
            ctx.moveTo(x0b, brY + 5); ctx.lineTo(x0b, brY);
            ctx.lineTo(x1b, brY); ctx.lineTo(x1b, brY + 5);
            ctx.stroke();
            // group caption centred above the rail
            const label = name === 'sigma+'
              ? `σ⁺ 组 · ${live.length} 条· 上能级 m_J ${fmtMList(g.origins)}`
              : `σ⁻ 组 · ${live.length} 条 · 上能级 m_J ${fmtMList(g.origins)}`;
            ctx.font = '11.5px system-ui, sans-serif';
            const tw = ctx.measureText(label).width;
            const cxm = (x0b + x1b) / 2;
            ctx.fillStyle = C.bg;
            ctx.fillRect(cxm - tw / 2 - 5, brY - 15, tw + 10, 15);
            ctx.fillStyle = colr; ctx.textAlign = 'center';
            ctx.fillText(label, cxm, brY - 4);
            ctx.textAlign = 'left';
            ctx.restore();
          });

        // a line tying the two groups together, showing they share one manifold
        ctx.save();
        ctx.font = '11px system-ui, sans-serif';
        const note = '两组均来自 6s7s ³S₁ 的三个 m_J 子能级；π 线沿 B 方向不辐射';
        ctx.fillStyle = C.dim; ctx.textAlign = 'center';
        ctx.fillText(note, ML + pw / 2, h - 30);
        ctx.textAlign = 'left';
        ctx.restore();
      }
    }

    // B = 0 annotation: single line
    if (B < 1e-6) {
      ctx.fillStyle = C.bg; ctx.globalAlpha = 0.55;
      ctx.fillRect(ML, MT, pw, ph);
      ctx.globalAlpha = 1;
      ctx.fillStyle = C.text; ctx.font = '15px system-ui, sans-serif';
      ctx.textAlign = 'center';
      ctx.fillText('B = 0：无场，546.1 nm 单条谱线', ML + pw / 2, MT + ph / 2);
      ctx.textAlign = 'left';
    }
    return { comps, smax };
  }

  /* =======================================================================
     3. Fabry-Perot interference rings.
        For each wavelength we accumulate the Airy transmission over angle,
        giving the real ring pattern that a CCD would record.
     ======================================================================= */
  function rings(canvas, B, opts) {
    opts = opts || {};
    const { ctx, w, h, dpr } = setupCanvas(canvas);
    const comps = P.applyGeometry(P.components(B), opts.view || 'perp', P.normalisePolariser(opts.polariser), P.normaliseQWP(opts.qwp));
    const fp = P.FP;

    ctx.fillStyle = C.bg; ctx.fillRect(0, 0, w, h);
    const cx = w / 2, cy = h / 2;
    // The panel is wider than it is tall, so a square field of view would leave
    // empty margins at the sides.  Map the full box to the field of view with a
    // single isotropic px-per-mm scale, so the fringes stay circular while the
    // image fills the available area.
    const thetaMaxDeg = opts.thetaMaxDeg || 4.0;
    const thMax = thetaMaxDeg * Math.PI / 180;
    const rMaxMm = fp.f_mm * Math.tan(thMax);
    const R = Math.min(w, h) / 2 - 6;
    const scalePxPerMm = Math.min((w / 2 - 4) / rMaxMm, (h / 2 - 4) / rMaxMm);

    const d_cm = fp.d_mm / 10;
    const active = comps.filter(c => c.weight >= P.VISIBLE_EPS);
    const norm = active.reduce((a, c) => a + c.strength * c.weight, 0) || 1;

    /* The pattern depends only on radius, so evaluate it once on a fine radial
       grid and paint each pixel by lookup.  Same mapping as radialProfile(),
       so what is displayed is exactly what gets measured. */
    const NR = 2048;
    const prof = new Float64Array(NR);
    for (let i = 0; i < NR; i++) {
      const rMm = (i / (NR - 1)) * rMaxMm;
      const ct = Math.cos(Math.atan(rMm / fp.f_mm));
      let v = 0;
      for (let k = 0; k < active.length; k++) {
        const c = active[k];
        const lam = c.lambda_nm * 1e-7;          // cm
        // half round-trip phase; see physics.airy() for the factor 2
        const half = Math.PI * 2 * fp.n * d_cm * ct / lam;
        const s = Math.sin(half);
        v += c.strength * c.weight / (1 + fp.F * s * s);
      }
      prof[i] = v / norm;
    }

    // Build the bitmap at BACKING resolution: putImageData ignores the ctx
    // transform, so an ImageData sized in CSS pixels would land in the top-left
    // corner of a dpr-scaled canvas and leave the rest blank.
    const bw = canvas.width, bh = canvas.height;
    const img = ctx.createImageData(bw, bh);
    const data = img.data;
    const cxB = bw / 2, cyB = bh / 2, RB = R * dpr, scaleB = scalePxPerMm * dpr;

    for (let py = 0; py < bh; py++) {
      const dy = py - cyB;
      for (let px = 0; px < bw; px++) {
        const dx = px - cxB;
        const rPix = Math.sqrt(dx * dx + dy * dy);
        let g = 0;
        if (rPix < RB) {
          const t = rPix / scaleB / rMaxMm;          // 0..1 across the field
          let idxf = t * (NR - 1);
          if (idxf < 0) idxf = 0; else if (idxf > NR - 1) idxf = NR - 1;
          const i0 = Math.floor(idxf), fr = idxf - i0;
          const val = prof[i0] * (1 - fr) + prof[Math.min(NR - 1, i0 + 1)] * fr;
          g = Math.pow(Math.min(1, val), 0.42);
          const rr = rPix / RB;
          g *= (1 - 0.55 * rr * rr);                  // radial vignette
        }
        const idx = (py * bw + px) * 4;
        data[idx] = Math.round(255 * Math.min(1, g * 0.92));
        data[idx + 1] = Math.round(255 * Math.min(1, g));
        data[idx + 2] = Math.round(255 * Math.min(1, 0.55 + 0.45 * g));
        data[idx + 3] = 255;
      }
    }
    ctx.save();
    ctx.setTransform(1, 0, 0, 1, 0, 0);
    ctx.putImageData(img, 0, 0);
    ctx.restore();

    // scale bar + labels (uses the same px-per-mm scale as the image)
    const barMm = opts.barMm || 2;
    const barPx = barMm * scalePxPerMm;
    ctx.strokeStyle = 'rgba(255,255,255,0.75)'; ctx.lineWidth = 2;
    ctx.beginPath();
    ctx.moveTo(14, h - 16); ctx.lineTo(14 + barPx, h - 16);
    ctx.moveTo(14, h - 21); ctx.lineTo(14, h - 11);
    ctx.moveTo(14 + barPx, h - 21); ctx.lineTo(14 + barPx, h - 11);
    ctx.stroke();
    ctx.fillStyle = '#fff'; ctx.font = '12px ui-monospace, monospace';
    ctx.fillText(`${barMm} mm`, 14 + barPx + 8, h - 12);

    // field-of-view note
    ctx.fillStyle = 'rgba(220,230,255,0.75)'; ctx.font = '12px system-ui, sans-serif';
    ctx.fillText(`视场 ±${thetaMaxDeg}°　θ_max`, 14, 20);

    /* ---- line-index annotation for the measurement ------------------------
       The lab manual numbers the nine rings of one interference order 1..9
       from the INNERMOST ring outwards, so a student can be told "measure
       rings 4 and 6" and actually find them.  When the caller supplies the
       order being measured (opts.markOrder), the nine radii are labelled
       along a horizontal ray, and the selected symmetric pair is picked out.
       Offsets are in the same px-per-mm scale as the image, so a label always
       sits on its own ring.                                                  */
    if (opts.markOrder && B > 1e-6) {
      const m = opts.markOrder;
      const pair = opts.markPair;                 // {inner, outer, d0}
      const lab = [];
      for (let i = 1; i <= 9; i++) {
        const D = P.splitDiameter(m, i, B, fp);
        if (!isFinite(D)) continue;
        lab.push({ i, r: D / 2, delta: P.deltaOfLineIndex(i) });
      }
      lab.sort((a, b) => a.r - b.r);              // innermost first
      if (lab.length) {
        /* Lay the labels along a horizontal ray to the right.  Successive
           rings of one order are only tens of micrometres apart, so the
           labels are stacked vertically with a guaranteed minimum gap and a
           leader line back to the exact radius — otherwise they overlap into
           an unreadable smear. */
        const ang = 0;                             // straight right
        const dx = Math.cos(ang), dy = Math.sin(ang);
        const rPx = lab.map(p => p.r * scalePxPerMm);
        const rMaxPx = rPx[rPx.length - 1];
        ctx.save();
        ctx.font = '12px ui-monospace, Menlo, monospace';
        ctx.textAlign = 'left';
        ctx.textBaseline = 'middle';
        /* Park the stack in the empty margin beside the pattern.  Only the
           nine rings of the measured order get labels, so the leader lines
           converge from a narrow radial band; they are drawn faintly except
           for the selected pair so the pattern stays readable. */
        const stackX = Math.max(cx + rMaxPx + 34, w * 0.5 + 40);
        const colW = 78;
        const x0 = Math.min(stackX, w - colW - 6);
        const span = 9 * 16;
        let y0 = cy - span / 2;
        y0 = Math.max(y0, 30);
        if (y0 + span > h - 8) y0 = Math.max(30, h - 8 - span);
        // panel behind the stack so it stays readable over any background
        ctx.fillStyle = 'rgba(8,12,24,0.62)';
        ctx.fillRect(x0 - 22, y0 - 22, colW + 28, span + 30);
        lab.forEach((p, n) => {
          const px = cx + rPx[n] * dx, py = cy + rPx[n] * dy;
          const ly = y0 + n * 16;
          const hot = pair && (p.i === pair.inner || p.i === pair.outer);
          const col = hot ? '#fbbf24' : 'rgba(235,242,255,0.90)';
          // leader from the ring to its label row
          ctx.strokeStyle = hot ? 'rgba(251,191,36,0.45)' : 'rgba(235,242,255,0.16)';
          ctx.lineWidth = hot ? 1.4 : 1;
          ctx.beginPath();
          ctx.moveTo(px + 5, py);
          ctx.lineTo(x0 - 20, ly);
          ctx.stroke();
          // dot on the ring itself
          ctx.fillStyle = col;
          ctx.beginPath(); ctx.arc(px, py, hot ? 3.6 : 2.2, 0, 6.284); ctx.fill();
          // index chip
          ctx.fillStyle = hot ? 'rgba(251,191,36,0.95)' : 'rgba(16,24,44,0.92)';
          ctx.beginPath(); ctx.arc(x0 - 11, ly, 9.5, 0, 6.284); ctx.fill();
          ctx.fillStyle = hot ? '#0a1020' : col;
          ctx.textAlign = 'center';
          ctx.fillText(String(p.i), x0 - 11, ly + 1);
          ctx.textAlign = 'left';
          // delta value
          ctx.fillStyle = hot ? 'rgba(251,191,36,0.85)' : 'rgba(200,214,240,0.66)';
          ctx.font = '10.5px ui-monospace, Menlo, monospace';
          ctx.fillText((p.delta > 0 ? '+' : '') + p.delta.toFixed(1), x0 + 1, ly + 1);
          ctx.font = '12px ui-monospace, Menlo, monospace';
        });
        // group caption
        ctx.font = '11.5px system-ui, sans-serif';
        const cap = '第 1…9 条 · m = ' + m;
        ctx.fillStyle = 'rgba(8,12,24,0.90)';
        ctx.fillRect(x0 - 22, y0 - 20, colW + 20, 16);
        ctx.fillStyle = 'rgba(235,242,255,0.92)';
        ctx.textAlign = 'left';
        ctx.fillText(cap, x0 - 19, y0 - 8);
        ctx.textBaseline = 'alphabetic';
        ctx.restore();
      }
    }

    if (B < 1e-6) {
      ctx.fillStyle = 'rgba(11,16,32,0.72)';
      const bw = 320, bh = 40, bx = (w - bw) / 2, by = h - 62;
      ctx.fillRect(bx, by, bw, bh);
      ctx.strokeStyle = C.dim; ctx.strokeRect(bx, by, bw, bh);
      ctx.fillStyle = C.text; ctx.font = '14px system-ui, sans-serif';
      ctx.textAlign = 'center';
      ctx.fillText('B = 0：无分裂，各波长干涉环重合', w / 2, by + 25);
      ctx.textAlign = 'left';
    }
    return { comps: active };
  }

  /* =======================================================================
     4. Radial intensity profile — used by the measurement tool.
        Returns an array of {r_mm, I, peaks:[...] } samples.
     ======================================================================= */
  function radialProfile(B, opts) {
    opts = opts || {};
    const fp = P.FP;
    const comps = P.applyGeometry(P.components(B), opts.view || 'perp', P.normalisePolariser(opts.polariser), P.normaliseQWP(opts.qwp))
      .filter(c => c.weight >= P.VISIBLE_EPS);
    const thMax = (opts.thetaMaxDeg || 4.0) * Math.PI / 180;
    const rMax = fp.f_mm * Math.tan(thMax);
    const N = opts.N || 1400;
    const d_cm = fp.d_mm / 10;

    const norm = comps.reduce((a, c) => a + c.strength * c.weight, 0) || 1;
    const prof = new Float64Array(N);
    for (let i = 0; i < N; i++) {
      const rMm = (i / (N - 1)) * rMax;
      const ct = Math.cos(Math.atan(rMm / fp.f_mm));
      let v = 0;
      for (let k = 0; k < comps.length; k++) {
        const c = comps[k];
        const lam = c.lambda_nm * 1e-7;
        const half = Math.PI * 2 * fp.n * d_cm * ct / lam;
        const s = Math.sin(half);
        v += c.strength * c.weight / (1 + fp.F * s * s);
      }
      prof[i] = v / norm;
    }
    return { prof, rMax, N };
  }

  /**
   * Find bright rings in a radial profile.
   *
   * Two distinct issues have to be handled:
   *  1. a high-finesse Airy ring is very narrow, so a single ring can trigger
   *     several adjacent local maxima -> cluster those together;
   *  2. genuine neighbouring rings are separated by a few tenths of a mm,
   *     which is much larger than the within-ring spacing but SMALLER than the
   *     overall span, so any tolerance derived from the total span will merge
   *     real rings.
   *
   * We therefore split the maxima into clusters using the MEDIAN gap, which is
   * dominated by the (many) within-ring steps rather than by the few
   * between-ring steps.
   */
  function findRings(profile, opts) {
    opts = opts || {};
    const { prof, rMax, N } = profile;
    const minProm = opts.minProm || 0.02;

    // 1. local maxima with sub-sample refinement
    const raw = [];
    for (let i = 1; i < N - 1; i++) {
      const v = prof[i];
      if (v > prof[i - 1] && v >= prof[i + 1]) {
        const a = prof[i - 1], b = prof[i], c = prof[i + 1];
        const den = a - 2 * b + c;
        const shift = Math.abs(den) > 1e-18 ? 0.5 * (a - c) / den : 0;
        const idx = i + Math.max(-1, Math.min(1, shift));
        raw.push({ r: (idx / (N - 1)) * rMax, I: v });
      }
    }
    if (!raw.length) return [];
    raw.sort((p, q) => p.r - q.r);

    // 2. cluster.  Within-ring maxima sit a fraction of a mm apart; the
    //    sample spacing rMax/N is the natural lower bound for that.
    const sampleStep = rMax / (N - 1);
    // a single Airy ring cannot be narrower than a few samples
    const mergeTol = Math.max(sampleStep * 3, rMax * 0.002);
    const merged = [raw[0]];
    for (let i = 1; i < raw.length; i++) {
      const last = merged[merged.length - 1];
      if (raw[i].r - last.r < mergeTol) {
        if (raw[i].I > last.I) merged[merged.length - 1] = raw[i];
      } else {
        merged.push(raw[i]);
      }
    }

    // 3. reject weak features against the OFF-RING background (the median),
    //    not against the dimmest ring.
    const floorLevel = opts.floorLevel != null ? opts.floorLevel : background(profile);
    return merged.filter(p => p.I - floorLevel > minProm);
  }

  /** Robust estimate of the off-ring background: median of the profile. */
  function background(profile) {
    const s = Array.from(profile.prof).sort((a, b) => a - b);
    return s[Math.floor(s.length / 2)] || 0;
  }

  /* ----------------------------------------------------------------- format */
  function fmtSci(x, digits) {
    if (x === 0) return '0';
    const d = digits === undefined ? 3 : digits;
    return x.toExponential(d).replace('e', ' × 10^').replace('+', '');
  }
  function fmtSigned(x, d) {
    const s = x.toFixed(d);
    return (x >= 0 ? '+' : '') + s;
  }
  /** Format a list of m_J values, e.g. [-1,0,1] -> "−1 · 0 · +1". */
  function fmtMList(arr) {
    return arr.map(v => (v > 0 ? '+' : v < 0 ? '−' : '') + Math.abs(v)).join(' · ');
  }

  root.ZeemanRender = {
    C, polColor, setupCanvas, energySVG, spectrum, rings,
    radialProfile, findRings, background, fmtSci
  };
})(typeof self !== 'undefined' ? self : this);
