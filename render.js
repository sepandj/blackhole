import * as P from './physics.js';
import { VERT, TRACE, DOWN, UP, FINAL } from './shaders.js';

const canvas = document.getElementById('c');
const gl = canvas.getContext('webgl2', { antialias: false, alpha: false, depth: false, stencil: false, powerPreference: 'high-performance', preserveDrawingBuffer: false });
if (!gl) { document.body.classList.add('nogl'); throw new Error('webgl2 not available'); }

const floatOK = !!gl.getExtension('EXT_color_buffer_float');
const halfOK = floatOK || !!gl.getExtension('EXT_color_buffer_half_float');
const FMT = halfOK ? { internal: gl.RGBA16F, type: gl.HALF_FLOAT } : { internal: gl.RGBA8, type: gl.UNSIGNED_BYTE };

// ---- programs ----
function compile(type, src) {
  const s = gl.createShader(type); gl.shaderSource(s, src); gl.compileShader(s);
  if (!gl.getShaderParameter(s, gl.COMPILE_STATUS)) throw new Error(gl.getShaderInfoLog(s));
  return s;
}
function program(fs) {
  const p = gl.createProgram();
  gl.attachShader(p, compile(gl.VERTEX_SHADER, VERT));
  gl.attachShader(p, compile(gl.FRAGMENT_SHADER, fs));
  gl.bindAttribLocation(p, 0, 'aPos');
  gl.linkProgram(p);
  if (!gl.getProgramParameter(p, gl.LINK_STATUS)) throw new Error(gl.getProgramInfoLog(p));
  const u = {};
  const n = gl.getProgramParameter(p, gl.ACTIVE_UNIFORMS);
  for (let i = 0; i < n; i++) { const a = gl.getActiveUniform(p, i); u[a.name] = gl.getUniformLocation(p, a.name); }
  return { p, u };
}
const progTrace = program(TRACE), progDown = program(DOWN), progUp = program(UP), progFinal = program(FINAL);

const vbo = gl.createBuffer();
gl.bindBuffer(gl.ARRAY_BUFFER, vbo);
gl.bufferData(gl.ARRAY_BUFFER, new Float32Array([-1, -1, 3, -1, -1, 3]), gl.STATIC_DRAW);
gl.enableVertexAttribArray(0);
gl.vertexAttribPointer(0, 2, gl.FLOAT, false, 0, 0);

// ---- lookup tables computed from the physics ----
function dataTex(w, internal, format, data) {
  const t = gl.createTexture();
  gl.bindTexture(gl.TEXTURE_2D, t);
  gl.texImage2D(gl.TEXTURE_2D, 0, internal, w, 1, 0, format, gl.FLOAT, data);
  gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.NEAREST);
  gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.NEAREST);
  gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
  gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
  return t;
}
const diskTable = P.diskTable();
const bbTable = P.blackbodyTable();
// float textures with linear filtering are not everywhere; the tables are fine enough for nearest
const texDisk = dataTex(P.DISK_N, gl.R32F, gl.RED, diskTable);
const texBB = dataTex(P.BB_N, gl.RGBA32F, gl.RGBA, bbTable);
const peak = P.diskPeak();
const tauPeak = Math.pow(peak.f, 0.25);
// luminance of a blackbody at temperature T, same table as the shader
function bbLum(T) {
  const u = Math.min(1, Math.max(0, Math.log(T / P.BB_TMIN) / Math.log(P.BB_TMAX / P.BB_TMIN)));
  const x = u * (P.BB_N - 1), i = Math.min(Math.floor(x), P.BB_N - 2), k = x - i;
  const at = (j) => bbTable[i * 4 + j] * (1 - k) + bbTable[(i + 1) * 4 + j] * k;
  return { rgb: [at(0), at(1), at(2)], lum: Math.pow(10, at(3)) };
}

// ---- render targets ----
function target(w, h) {
  const t = gl.createTexture();
  gl.bindTexture(gl.TEXTURE_2D, t);
  gl.texImage2D(gl.TEXTURE_2D, 0, FMT.internal, w, h, 0, gl.RGBA, FMT.type, null);
  gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.LINEAR);
  gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.LINEAR);
  gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
  gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
  const f = gl.createFramebuffer();
  gl.bindFramebuffer(gl.FRAMEBUFFER, f);
  gl.framebufferTexture2D(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT0, gl.TEXTURE_2D, t, 0);
  return { t, f, w, h };
}
let levels = [], tw = 0, th = 0;
function free() { for (const l of levels) { gl.deleteTexture(l.t); gl.deleteFramebuffer(l.f); } levels = []; }
function build(w, h) {
  free();
  tw = w; th = h;
  levels.push(target(w, h));
  let lw = w, lh = h;
  for (let i = 0; i < 6 && lw > 8 && lh > 8; i++) { lw = Math.max(1, lw >> 1); lh = Math.max(1, lh >> 1); levels.push(target(lw, lh)); }
}

// ---- view state ----
const coarse = matchMedia('(pointer: coarse)').matches;
const params = new URLSearchParams(location.search);
const state = {
  hole: P.HOLES[0],
  edd: 2e-4,
  yaw: 0.6, pitch: 0.12, dist: 70, target: 22,
  vyaw: 0, vpitch: 0,
  t: 0,
  expo: 0, // adapts
  lastInput: -1e9,
  intro: true,
};
if (params.has('d')) { state.dist = state.target = +params.get('d'); state.intro = false; }
if (params.has('p')) state.pitch = +params.get('p') * Math.PI / 180;
if (params.has('y')) state.yaw = +params.get('y') * Math.PI / 180;
if (params.has('h')) state.hole = P.HOLES.find((h) => h.id === params.get('h')) || state.hole;
if (params.has('e')) state.edd = Math.pow(10, +params.get('e'));
const STILL = params.has('still');
if (STILL) state.intro = false;
if (params.has('t')) state.t = +params.get('t');
const DMIN = 1.25, DMAX = 120;
const ISCO_PERIOD = 2 * Math.PI / P.orbitOmega(P.ISCO); // r_s / c
const SPEED = ISCO_PERIOD / 18; // coordinate time per real second, far away: one inner orbit in 18 s

// ---- input: drag to fly around, pinch or scroll to move in and out ----
const pointers = new Map();
let pinch0 = 0, dist0 = 0;
function onDown(e) {
  if (e.target !== canvas) return;
  canvas.setPointerCapture(e.pointerId);
  pointers.set(e.pointerId, { x: e.clientX, y: e.clientY, t: performance.now(), x0: e.clientX, y0: e.clientY, t0: performance.now(), moved: 0 });
  if (pointers.size > 1) for (const q of pointers.values()) q.moved = 1e9;
  state.vyaw = state.vpitch = 0;
  state.intro = false; state.target = state.dist;
  if (pointers.size === 2) { const [a, b] = [...pointers.values()]; pinch0 = Math.hypot(a.x - b.x, a.y - b.y); dist0 = state.dist; }
  state.lastInput = performance.now();
  document.body.classList.add('touched');
}
function onMove(e) {
  const p = pointers.get(e.pointerId); if (!p) return;
  const now = performance.now();
  const dx = e.clientX - p.x, dy = e.clientY - p.y;
  const dtp = Math.max(1, now - p.t);
  p.moved += Math.hypot(dx, dy);
  p.x = e.clientX; p.y = e.clientY; p.t = now;
  if (pointers.size === 1) {
    const k = 3.2 / Math.max(innerWidth, innerHeight) * Math.min(1, Math.sqrt(state.dist / 10) * 0.8 + 0.2);
    state.yaw -= dx * k; state.pitch += dy * k;
    state.vyaw = -dx * k / dtp * 1000; state.vpitch = dy * k / dtp * 1000;
  } else if (pointers.size === 2) {
    const [a, b] = [...pointers.values()];
    const d = Math.hypot(a.x - b.x, a.y - b.y);
    if (pinch0 > 0) state.dist = state.target = clampD(dist0 * pinch0 / d);
  }
  state.lastInput = now;
}
function onUp(e) {
  const p0 = pointers.get(e.pointerId);
  if (p0 && pointers.size === 1 && e.type === 'pointerup' && p0.moved < 9 && performance.now() - p0.t0 < 400) throwAt(e.clientX, e.clientY);
  pointers.delete(e.pointerId);
  if (pointers.size < 2) pinch0 = 0;
  if (pointers.size === 1) { const p = [...pointers.values()][0]; p.t = performance.now(); state.vyaw = state.vpitch = 0; }
}
const clampD = (d) => Math.min(DMAX, Math.max(DMIN, d));
canvas.addEventListener('pointerdown', onDown);
canvas.addEventListener('pointermove', onMove);
canvas.addEventListener('pointerup', onUp);
canvas.addEventListener('pointercancel', onUp);
canvas.addEventListener('wheel', (e) => {
  e.preventDefault();
  state.intro = false;
  state.target = clampD(state.target * Math.exp(e.deltaY * (e.deltaMode ? 0.05 : 0.0012)));
  state.lastInput = performance.now();
  document.body.classList.add('touched');
}, { passive: false });
addEventListener('keydown', (e) => {
  const k = { ArrowLeft: [0.08, 0], ArrowRight: [-0.08, 0], ArrowUp: [0, 0.06], ArrowDown: [0, -0.06] }[e.key];
  if (k) { state.yaw += k[0]; state.pitch += k[1]; state.lastInput = performance.now(); }
  if (e.key === '+' || e.key === '=') state.target = clampD(state.target / 1.15);
  if (e.key === '-') state.target = clampD(state.target * 1.15);
});

// ---- probes: tap to throw one ----
// thrown from where you hover, towards where you tapped, at the speed of a circular orbit
// at your distance. sideways it orbits; towards the hole it falls.
const PROBE_T = [6500, 3600, 15000];        // colour of each probe's lamp, kelvin
const probes = [];
let probeSlot = 0, view = null, camKey = 0, message = '', messageUntil = 0;
function throwAt(cx, cy) {
  if (!view) return;
  // how far from the hole you tap sets how sideways it goes: near the hole straight in,
  // near the edge of the screen (or further) sideways, which at this speed is a circle
  const W = canvas.clientWidth, H = canvas.clientHeight;
  const ox = cx - W / 2, oy = H / 2 - cy;
  const s = Math.min(1, Math.hypot(ox, oy) / (0.42 * Math.min(W, H)));
  const th = s * Math.PI / 2;
  const B = view.B, ol = Math.hypot(ox, oy) || 1;
  const side = [0, 1, 2].map((i) => (B[i] * ox + B[3 + i] * oy) / ol);
  const fw = [B[6], B[7], B[8]];
  const d = [0, 1, 2].map((i) => Math.cos(th) * fw[i] + Math.sin(th) * side[i]);
  // start a little inside your own radius so a full circle doesn't come back through you
  const x0 = view.cam.map((c) => c * 0.88);
  const r = Math.hypot(...x0);
  const w = r > 1.52 ? Math.min(0.97, P.circularSpeed(r)) : 0.97;
  const k = probeSlot++ % 3;
  const old = probes.findIndex((q) => q.k === k); if (old >= 0) probes.splice(old, 1);
  const T = PROBE_T[k];
  const c = bbLum(T), m = Math.max(...c.rgb);
  probes.push({ k, T, K: 1.5 / c.lum, rgb: c.rgb.map((v) => Math.round(255 * Math.min(1, v / m))), s: P.throwProbe(x0, d, w), w0: w, trail: [], gone: false });
  state.lastInput = performance.now();
  document.body.classList.add('touched', 'thrown');
}
function stepProbes(dtc) {
  for (let i = probes.length - 1; i >= 0; i--) {
    const q = probes[i];
    P.stepProbe(q.s, dtc, { k: 0.01, maxSub: 3000 });
    const r = Math.hypot(...q.s.x);
    const last = q.trail[q.trail.length - 1];
    if (!last || Math.hypot(q.s.x[0] - last.x[0], q.s.x[1] - last.x[1], q.s.x[2] - last.x[2]) > Math.max(0.06, 0.03 * r)) {
      q.trail.push({ x: q.s.x.slice(), psi: last ? last.psi : null, d: null, hidden: false, key: -1 });
      if (q.trail.length > 240) q.trail.shift();
    }
    const loc = P.probeLocal(q.s);
    if (loc.clock < 0.004 || r < 1.00002) {
      probes.splice(i, 1);
      message = 'the probe is gone. from out here it never quite crosses the horizon, it slows down, reddens and fades';
      messageUntil = performance.now() + 9000;
    } else if (r > 160) {
      probes.splice(i, 1);
      message = 'the probe got away';
      messageUntil = performance.now() + 5000;
    }
  }
}
// solve where each trail point's light lands, a few milliseconds per frame
let solveAt = 0;
function solveTrails(budgetMs) {
  const t0 = performance.now();
  const all = [];
  for (const q of probes) for (const pt of q.trail) all.push(pt);
  if (!all.length) return;
  for (let n = 0; n < all.length && performance.now() - t0 < budgetMs; n++) {
    const pt = all[(solveAt + n) % all.length];
    if (pt.key === camKey) continue;
    const im = P.imageOf(view.cam, pt.x, pt.psi);
    pt.key = camKey;
    if (im) { pt.d = im.d; pt.psi = im.psi; pt.hidden = im.hidden; } else pt.d = null;
  }
  solveAt = (solveAt + 7) % Math.max(1, all.length);
}
const ov = document.getElementById('o');
const octx = ov.getContext('2d');
function drawTrails() {
  const dpr = Math.min(devicePixelRatio || 1, 2);
  const W = Math.round(ov.clientWidth * dpr), H = Math.round(ov.clientHeight * dpr);
  if (ov.width !== W || ov.height !== H) { ov.width = W; ov.height = H; }
  octx.clearRect(0, 0, W, H);
  if (!probes.length || !view) return;
  const B = view.B;
  const proj = (d) => {
    const z = d[0] * B[6] + d[1] * B[7] + d[2] * B[8];
    if (z <= 0.01) return null;
    const x = (d[0] * B[0] + d[1] * B[1] + d[2] * B[2]) / z / view.tan[0];
    const y = (d[0] * B[3] + d[1] * B[4] + d[2] * B[5]) / z / view.tan[1];
    return [(x + 1) / 2 * W, (1 - y) / 2 * H];
  };
  octx.globalCompositeOperation = 'lighter';
  octx.lineCap = 'round';
  octx.lineWidth = 1.5 * dpr;
  const jump = 0.2 * Math.min(W, H);
  for (const q of probes) {
    const n = q.trail.length;
    let prev = null;
    for (let i = 0; i < n; i++) {
      const pt = q.trail[i];
      const sp = pt.d ? proj(pt.d) : null;
      if (sp && prev && Math.hypot(sp[0] - prev[0], sp[1] - prev[1]) < jump) {
        const a = 0.85 * Math.pow((i + 1) / n, 1.1) * (pt.hidden ? 0.25 : 1);
        octx.strokeStyle = `rgba(${q.rgb[0]},${q.rgb[1]},${q.rgb[2]},${a.toFixed(3)})`;
        octx.beginPath(); octx.moveTo(prev[0], prev[1]); octx.lineTo(sp[0], sp[1]); octx.stroke();
      }
      prev = sp;
    }
  }
}

// ---- resolution that adapts to the device ----
let scale = STILL ? 1 : coarse ? 0.55 : 0.8;
const dprCap = coarse ? 2 : 2;
let steps = coarse ? 260 : 360;
let frameTimes = [], lastAdjust = 0;
function adapt(now, dt) {
  frameTimes.push(dt);
  if (frameTimes.length < 24 || now - lastAdjust < 700) return;
  const s = frameTimes.slice().sort((a, b) => a - b);
  const med = s[s.length >> 1];
  frameTimes = [];
  lastAdjust = now;
  if (med > 26 && scale > 0.3) scale = Math.max(0.3, scale * 0.85);
  else if (med < 18.5 && scale < 1) scale = Math.min(1, scale * 1.06);
}

function size() {
  const dpr = Math.min(devicePixelRatio || 1, dprCap);
  const cw = Math.round(canvas.clientWidth * dpr), ch = Math.round(canvas.clientHeight * dpr);
  if (canvas.width !== cw || canvas.height !== ch) { canvas.width = cw; canvas.height = ch; }
  const w = Math.max(64, Math.round(cw * scale)), h = Math.max(64, Math.round(ch * scale));
  if (w !== tw || h !== th) build(w, h);
}

// ---- hud ----
const $ = (id) => document.getElementById(id);
const fmt = (x, d = 3) => {
  if (!isFinite(x)) return '-';
  if (Math.abs(x) >= 1e5 || (Math.abs(x) < 1e-2 && x !== 0)) {
    const e = Math.floor(Math.log10(Math.abs(x)));
    return `${(x / 10 ** e).toFixed(1)}×10${sup(e)}`;
  }
  return Number(x.toPrecision(d)).toLocaleString('en');
};
const sup = (n) => String(n).split('').map((c) => '⁰¹²³⁴⁵⁶⁷⁸⁹⁻'['0123456789-'.indexOf(c)]).join('');
function distText(m) {
  if (m < 1e9) return `${fmt(m / 1e3)} km`;
  if (m < 0.1 * P.LY) return `${fmt(m / P.AU)} au`;
  return `${fmt(m / P.LY)} light years`;
}
function durText(s) {
  if (s < 1e-3) return `${fmt(s * 1e6)} microseconds`;
  if (s < 1) return `${fmt(s * 1e3)} milliseconds`;
  if (s < 120) return `${fmt(s)} seconds`;
  if (s < 7200) return `${fmt(s / 60)} minutes`;
  if (s < 3 * 86400) return `${fmt(s / 3600)} hours`;
  if (s < 2 * P.YEAR) return `${fmt(s / 86400)} days`;
  return `${fmt(s / P.YEAR)} years`;
}
function massText(m) {
  if (m >= 1e9) return `${fmt(m / 1e9, 2)} billion suns`;
  if (m >= 1e6) return `${fmt(m / 1e6, 3)} million suns`;
  return `${fmt(m, 3)} suns`;
}
function hud() {
  const h = state.hole, m = h.mass * P.MSUN;
  const rs = P.schwarzschildRadius(m), tu = P.timeUnit(m);
  const T0 = P.diskT0(m, state.edd);
  const r = state.dist, f = P.clockRate(r);
  $('name').textContent = h.name;
  $('mass').textContent = massText(h.mass);
  $('l1').textContent = `${fmt(r, 3)} rₛ from the centre, ${distText(r * rs)}`;
  $('l2').textContent = r < 1.5
    ? `your clock runs at ${fmt(f * 100, 3)}% of a far one. light can't orbit this close`
    : `your clock runs at ${fmt(f * 100, 4)}% of a far one`;
  $('l3').textContent = `inner orbit ${durText(ISCO_PERIOD * tu)}, hottest gas ${fmt(Math.round(T0 * tauPeak / 100) * 100, 3)} K`;
  $('l4').textContent = `sped up, each second is ${durText(SPEED * warp * tu)} of your time`;
  const q = probes[probes.length - 1];
  if (q) {
    const loc = P.probeLocal(q.s);
    $('l5').textContent = `probe ${fmt(loc.r, 3)} r\u209b out, ${fmt(loc.speed, 2)} c, its clock at ${fmt(loc.clock * 100, 3)}%`;
  } else $('l5').textContent = performance.now() < messageUntil ? message : '';
  $('feed').textContent = state.edd >= 0.01 ? `${fmt(state.edd * 100, 2)}%` : `${fmt(state.edd * 100, 2)}%`;
}
$('name').addEventListener('click', () => {
  const i = P.HOLES.indexOf(state.hole);
  state.hole = P.HOLES[(i + 1) % P.HOLES.length];
  hud();
});
$('rate').addEventListener('input', (e) => { state.edd = Math.pow(10, +e.target.value); hud(); });
$('rate').value = Math.log10(state.edd);

// ---- frame ----
let last = performance.now(), fade = 0, warp = 1;
function basis(yaw, pitch) {
  const cp = Math.cos(pitch);
  const pos = [cp * Math.sin(yaw), Math.sin(pitch), cp * Math.cos(yaw)];
  const fw = pos.map((x) => -x);
  let rt = [fw[2], 0, -fw[0]]; // cross(fw, up=(0,1,0)) up to sign
  const rl = Math.hypot(rt[0], rt[2]) || 1; rt = rt.map((x) => x / rl);
  // make right-handed: right = fw x up
  rt = [-rt[0], 0, -rt[2]];
  const up = [rt[1] * fw[2] - rt[2] * fw[1], rt[2] * fw[0] - rt[0] * fw[2], rt[0] * fw[1] - rt[1] * fw[0]];
  return { pos, m: [...rt, ...up, ...fw] };
}
function frame(now) {
  const dt = Math.min(0.1, (now - last) / 1000);
  if (!STILL) adapt(now, (now - last));
  last = now;

  // motion
  if (pointers.size === 0) {
    state.yaw += state.vyaw * dt; state.pitch += state.vpitch * dt;
    const damp = Math.exp(-dt * 2.2);
    state.vyaw *= damp; state.vpitch *= damp;
    if (now - state.lastInput > 6000) state.yaw += dt * 0.018 * (1 - Math.exp(-(now - state.lastInput - 6000) / 4000));
  }
  state.pitch = Math.max(-1.553, Math.min(1.553, state.pitch));
  if (state.intro) state.target = 22;
  const ease = state.intro ? 0.35 : 4;
  state.dist += (state.target - state.dist) * (1 - Math.exp(-dt * ease));
  if (state.intro && Math.abs(state.dist - 22) < 0.05) state.intro = false;

  // physics for this frame
  const r0 = state.dist, f = Math.sqrt(1 - 1 / r0);
  // with a probe out, time runs faster so its orbit is watchable (one kepler period in ~16 s)
  const q = probes[probes.length - 1];
  const wantWarp = q ? Math.min(30, Math.max(1, 2 * Math.PI * Math.sqrt(Math.hypot(...q.s.x) ** 3 / 0.5) / 16 / SPEED)) : 1;
  warp += (wantWarp - warp) * (1 - Math.exp(-dt * 1.2));
  const dtc = STILL ? 0 : dt * SPEED * warp / f; // the hovering camera ages slower, so the gas turns faster to it
  state.t += dtc;
  if (dtc > 0) stepProbes(dtc);
  const m = state.hole.mass * P.MSUN;
  const T0 = P.diskT0(m, state.edd);
  // eyes adapt: aim the hottest gas, as seen from the hovering camera, at a fixed display level
  const want = 0.6 / bbLum(T0 * tauPeak / f).lum;
  state.expo = state.expo ? Math.exp(Math.log(state.expo) + (Math.log(want) - Math.log(state.expo)) * (1 - Math.exp(-dt * 2.5))) : want;
  const s0 = bbLum(5000), s1 = bbLum(5000 / f);
  const skyGain = [0, 1, 2].map((i) => (s1.rgb[i] * s1.lum) / (s0.rgb[i] * s0.lum));

  size();
  const { pos, m: B } = basis(state.yaw, state.pitch);
  const cam = pos.map((x) => x * r0);
  if (!view || Math.hypot(cam[0] - view.cam[0], cam[1] - view.cam[1], cam[2] - view.cam[2]) > 1e-4 * r0) camKey++;
  const aspect = tw / th;
  // look wider as the shadow grows (synge: sin a = (sqrt 27 / 2) / r * sqrt(1 - 1/r))
  const sa = Math.min(1, P.SHADOW_B / r0 * f);
  const shadowHalf = r0 > 1.5 ? Math.asin(sa) : Math.PI - Math.asin(sa);
  const shortHalf = Math.tan(Math.min(68, Math.max(aspect >= 1 ? 31 : 25, shadowHalf * 180 / Math.PI * 1.55)) * Math.PI / 180);
  const tan = aspect >= 1 ? [shortHalf * aspect, shortHalf] : [shortHalf, shortHalf / aspect];
  view = { cam, B, tan };

  gl.disable(gl.BLEND);
  gl.bindFramebuffer(gl.FRAMEBUFFER, levels[0].f);
  gl.viewport(0, 0, tw, th);
  const u = progTrace.u;
  gl.useProgram(progTrace.p);
  gl.uniform2f(u.uRes, tw, th);
  gl.uniform3fv(u.uCam, cam);
  gl.uniformMatrix3fv(u.uBasis, false, B);
  gl.uniform2fv(u.uTan, tan);
  gl.uniform1f(u.uF, f);
  gl.uniform1f(u.uTime, state.t % 1e5);
  gl.uniform1f(u.uT0, T0);
  gl.uniform1f(u.uExpo, state.expo);
  gl.uniform1f(u.uSky, 1.0);
  gl.uniform3fv(u.uSkyGain, skyGain);
  gl.uniform1f(u.uPix, 2 * tan[1] / th);
  gl.uniform1i(u.uSteps, steps);
  gl.uniform1f(u.uDiskOut, P.DISK_OUT);
  const pr = new Float32Array(12), pw = new Float32Array(12), pt = new Float32Array(6);
  for (const q of probes) {
    const loc = P.probeLocal(q.s), k = q.k;
    pr.set([...q.s.x, 1], k * 4); pw.set([...loc.w, loc.gamma], k * 4); pt.set([q.T, q.K], k * 2);
  }
  gl.uniform4fv(u['uPr[0]'], pr); gl.uniform4fv(u['uPw[0]'], pw); gl.uniform2fv(u['uPt[0]'], pt);
  gl.activeTexture(gl.TEXTURE0); gl.bindTexture(gl.TEXTURE_2D, texDisk); gl.uniform1i(u.uDiskT, 0);
  gl.activeTexture(gl.TEXTURE1); gl.bindTexture(gl.TEXTURE_2D, texBB); gl.uniform1i(u.uBB, 1);
  gl.drawArrays(gl.TRIANGLES, 0, 3);

  // bloom
  gl.useProgram(progDown.p);
  gl.activeTexture(gl.TEXTURE0);
  for (let i = 1; i < levels.length; i++) {
    const s = levels[i - 1], d = levels[i];
    gl.bindFramebuffer(gl.FRAMEBUFFER, d.f); gl.viewport(0, 0, d.w, d.h);
    gl.bindTexture(gl.TEXTURE_2D, s.t);
    gl.uniform1i(progDown.u.uSrc, 0);
    gl.uniform2f(progDown.u.uTexel, 1 / s.w, 1 / s.h);
    gl.drawArrays(gl.TRIANGLES, 0, 3);
  }
  gl.useProgram(progUp.p);
  gl.enable(gl.BLEND); gl.blendFunc(gl.ONE, gl.ONE);
  for (let i = levels.length - 1; i > 1; i--) {
    const s = levels[i], d = levels[i - 1];
    gl.bindFramebuffer(gl.FRAMEBUFFER, d.f); gl.viewport(0, 0, d.w, d.h);
    gl.bindTexture(gl.TEXTURE_2D, s.t);
    gl.uniform1i(progUp.u.uSrc, 0);
    gl.uniform2f(progUp.u.uTexel, 1 / s.w, 1 / s.h);
    gl.drawArrays(gl.TRIANGLES, 0, 3);
  }
  gl.disable(gl.BLEND);

  gl.bindFramebuffer(gl.FRAMEBUFFER, null);
  gl.viewport(0, 0, canvas.width, canvas.height);
  gl.useProgram(progFinal.p);
  gl.activeTexture(gl.TEXTURE0); gl.bindTexture(gl.TEXTURE_2D, levels[0].t); gl.uniform1i(progFinal.u.uScene, 0);
  gl.activeTexture(gl.TEXTURE1); gl.bindTexture(gl.TEXTURE_2D, levels[Math.min(1, levels.length - 1)].t); gl.uniform1i(progFinal.u.uBloom, 1);
  gl.uniform1f(progFinal.u.uBloomK, 0.14 / Math.max(1, levels.length - 2));
  fade = STILL ? 1 : Math.min(1, fade + dt / 2.5);
  gl.uniform1f(progFinal.u.uFade, fade * fade * (3 - 2 * fade));
  gl.uniform2f(progFinal.u.uRes, canvas.width, canvas.height);
  gl.drawArrays(gl.TRIANGLES, 0, 3);

  if (probes.length) solveTrails(coarse ? 2.5 : 3.5);
  drawTrails();
  if (now - (frame.h || 0) > 250) { frame.h = now; hud(); }
  if (STILL) { gl.finish(); window.__bh.done = true; return; }
  requestAnimationFrame(frame);
}
hud();
requestAnimationFrame(frame);

// for testing
window.__bh = { state, probes, throwAt, get warp() { return warp; }, get scale() { return scale; }, set scale(v) { scale = v; }, set steps(v) { steps = v; } };
