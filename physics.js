// physics for the black hole. SI units unless a name says otherwise.
// the renderer works in units of the schwarzschild radius (r_s = 1) and
// light-crossing times of it (r_s / c); everything here converts to and from those.

export const G = 6.67430e-11;        // m^3 kg^-1 s^-2
export const C = 299792458;          // m / s
export const MSUN = 1.98892e30;      // kg
export const SIGMA = 5.670374419e-8; // stefan-boltzmann, W m^-2 K^-4
export const H = 6.62607015e-34;     // planck, J s
export const KB = 1.380649e-23;      // boltzmann, J / K
export const MP = 1.67262192e-27;    // proton mass, kg
export const SIGMA_T = 6.6524587e-29;// thomson cross section, m^2
export const AU = 1.495978707e11;    // m
export const LY = 9.4607304725808e15;// m
export const YEAR = 365.25 * 86400;  // s

// three real black holes, masses from the literature.
export const HOLES = [
  { id: 'm87', name: 'm87*', mass: 6.5e9, where: 'centre of the galaxy messier 87, 55 million light years away', ref: 'event horizon telescope 2019' },
  { id: 'sgra', name: 'sagittarius a*', mass: 4.297e6, where: 'centre of our milky way, 27,000 light years away', ref: 'gravity collaboration 2019' },
  { id: 'cyg', name: 'cygnus x-1', mass: 21.2, where: 'a star-eating black hole in our galaxy, 7,200 light years away', ref: 'miller-jones et al. 2021' },
];

export const schwarzschildRadius = (mKg) => 2 * G * mKg / (C * C);   // m
export const timeUnit = (mKg) => schwarzschildRadius(mKg) / C;        // s, one r_s / c

// a static (hovering) observer at radius r (in r_s) ages slower than one far away
// by sqrt(1 - r_s / r).
export const clockRate = (r) => Math.sqrt(1 - 1 / r);

// circular orbits in schwarzschild, radius in r_s, time in r_s / c.
// omega = sqrt(G M / r^3) in coordinate time; with r_s = 1, G M = 1/2.
export const orbitOmega = (r) => Math.sqrt(0.5 / (r * r * r));
export const ISCO = 3;          // innermost stable circular orbit, r_s
export const PHOTON_SPHERE = 1.5;
export const SHADOW_B = Math.sqrt(27) / 2;   // critical impact parameter, r_s

// ---- thin disk (novikov-thorne / page-thorne, schwarzschild) ----
// flux from one face in geometric units (G = c = M = 1, mdot = 1), radius in M.
// F = mdot / (4 pi r) * (-dOmega/dr) / (E - Omega L)^2 * integral_isco^r (E - Omega L) dL/dr dr
function ptParts(r) {
  const s = Math.sqrt(1 - 3 / r);
  const E = (1 - 2 / r) / s;
  const L = Math.sqrt(r) / s;
  const Om = Math.pow(r, -1.5);
  return { E, L, Om };
}
export function pageThorneFluxM(rM) {
  if (rM <= 6) return 0;
  const n = 2000, a = 6, b = rM;
  let sum = 0;
  // simpson on (E - Omega L) dL/dr, dL/dr by central difference
  const f = (r) => {
    const { E, L, Om } = ptParts(r);
    const d = 1e-5 * r;
    const dL = (ptParts(r + d).L - ptParts(Math.max(6.0000001, r - d)).L) / (r + d - Math.max(6.0000001, r - d));
    return (E - Om * L) * dL;
  };
  const hh = (b - a) / n;
  for (let i = 0; i <= n; i++) {
    const r = a + i * hh;
    const w = i === 0 || i === n ? 1 : (i % 2 ? 4 : 2);
    sum += w * f(Math.max(r, 6.0000001));
  }
  const integral = sum * hh / 3;
  const { E, L, Om } = ptParts(rM);
  const dOm = 1.5 * Math.pow(rM, -2.5);
  return (1 / (4 * Math.PI * rM)) * dOm / ((E - Om * L) ** 2) * integral;
}

// accretion. mdot as a fraction of the eddington rate.
export const ETA = 1 - Math.sqrt(8 / 9);   // radiative efficiency of a schwarzschild thin disk
export const eddingtonLuminosity = (mKg) => 4 * Math.PI * G * mKg * MP * C / SIGMA_T; // W
export const eddingtonMdot = (mKg) => eddingtonLuminosity(mKg) / (ETA * C * C);      // kg / s

// temperature scale T0 such that T(r) = T0 * F_M(r / M)^(1/4)
export function diskT0(mKg, eddFrac) {
  const mdot = eddFrac * eddingtonMdot(mKg);
  const rg = G * mKg / (C * C);
  return Math.pow(mdot * C * C / (rg * rg * SIGMA), 0.25);
}

// table of F^(1/4) over the disk, radius in r_s, sampled in log r
export const DISK_IN = ISCO, DISK_OUT = 30, DISK_N = 256;
export function diskTable() {
  const t = new Float32Array(DISK_N);
  for (let i = 0; i < DISK_N; i++) {
    const r = DISK_IN * Math.pow(DISK_OUT / DISK_IN, i / (DISK_N - 1));
    t[i] = Math.pow(pageThorneFluxM(2 * r), 0.25);
  }
  return t;
}
export function diskPeak() {
  let best = 0, at = 0;
  for (let i = 0; i < 4000; i++) {
    const r = 3 + i * 0.005;
    const f = pageThorneFluxM(2 * r);
    if (f > best) { best = f; at = r; }
  }
  return { r: at, f: best };
}

// ---- blackbody colour as the eye would see it ----
// cie 1931 matching functions, multi-lobe fit (wyman, sloan, shirley 2013)
function lobe(l, mu, s1, s2) { const t = (l - mu) / (l < mu ? s1 : s2); return Math.exp(-0.5 * t * t); }
const xbar = (l) => 1.056 * lobe(l, 599.8, 37.9, 31.0) + 0.362 * lobe(l, 442.0, 16.0, 26.7) - 0.065 * lobe(l, 501.1, 20.4, 26.2);
const ybar = (l) => 0.821 * lobe(l, 568.8, 46.9, 40.5) + 0.286 * lobe(l, 530.9, 16.3, 31.1);
const zbar = (l) => 1.217 * lobe(l, 437.0, 11.8, 36.0) + 0.681 * lobe(l, 459.0, 26.0, 13.8);

// planck spectral radiance B_lambda, W sr^-1 m^-3
function planck(lm, T) {
  const x = H * C / (lm * KB * T);
  if (x > 700) return 0;
  return (2 * H * C * C) / (lm ** 5) / Math.expm1(x);
}

// returns linear srgb chromaticity (scaled so luminance = 1) and luminance in cd/m^2
export function blackbody(T) {
  let X = 0, Y = 0, Z = 0;
  for (let l = 360; l <= 830; l += 1) {
    const b = planck(l * 1e-9, T) * 1e-9; // per nm
    X += b * xbar(l); Y += b * ybar(l); Z += b * zbar(l);
  }
  const lum = 683 * Y;
  if (Y <= 0) return { rgb: [1, 0, 0], lum: 0 };
  X /= Y; Z /= Y; Y = 1;
  let r = 3.2406 * X - 1.5372 * Y - 0.4986 * Z;
  let g = -0.9689 * X + 1.8758 * Y + 0.0415 * Z;
  let b = 0.0557 * X - 0.2040 * Y + 1.0570 * Z;
  // out of gamut (deep red): pull towards grey at the same luminance
  const m = Math.min(r, g, b);
  if (m < 0) {
    const k = -m / (1 - m);
    r = r + (1 - r) * k; g = g + (1 - g) * k; b = b + (1 - b) * k;
    const y = 0.2126 * r + 0.7152 * g + 0.0722 * b;
    r /= y; g /= y; b /= y;
  }
  return { rgb: [r, g, b], lum };
}

export const BB_TMIN = 300, BB_TMAX = 1e9, BB_N = 512;
export function blackbodyTable() {
  const t = new Float32Array(BB_N * 4);
  for (let i = 0; i < BB_N; i++) {
    const T = BB_TMIN * Math.pow(BB_TMAX / BB_TMIN, i / (BB_N - 1));
    const { rgb, lum } = blackbody(T);
    t[i * 4] = rgb[0]; t[i * 4 + 1] = rgb[1]; t[i * 4 + 2] = rgb[2];
    t[i * 4 + 3] = Math.log10(Math.max(lum, 1e-30));
  }
  return t;
}

// ---- photon paths, the same integrator the shader runs ----
// in units r_s = 1 a light ray in schwarzschild obeys, exactly,
//   d2x/dl2 = -3/2 h^2 x / r^5,  h = |x cross dx/dl|
// (equivalent to u'' + u = 3/2 u^2 for u = 1/r along the orbit plane)
export function accel(x, h2) {
  const r2 = x[0] * x[0] + x[1] * x[1] + x[2] * x[2];
  const k = -1.5 * h2 / (r2 * r2 * Math.sqrt(r2));
  return [k * x[0], k * x[1], k * x[2]];
}
const add = (a, b, s) => [a[0] + b[0] * s, a[1] + b[1] * s, a[2] + b[2] * s];
const cross = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
const len = (a) => Math.sqrt(dot(a, a));

// remaining bend from a point in the weak field out to infinity, to first order in r_s / r
export function tailBend(x, v) {
  const r = len(x), vl = len(v);
  const cps = dot(x, v) / (r * vl);
  const b = len(cross(x, v)) / vl;
  return (1 - cps * (3 - cps * cps) / 2) / b;
}

// step control shared with the shader
export const STEP_K = 0.075, MAX_STEPS = 400;
export function trace(x0, v0, opts = {}) {
  const k = opts.k ?? STEP_K, maxSteps = opts.maxSteps ?? MAX_STEPS;
  let x = x0.slice(), v = v0.slice();
  const h2 = dot(cross(x, v), cross(x, v));
  const rEsc = Math.max(1.1 * len(x0), 40);
  for (let i = 0; i < maxSteps; i++) {
    const r = len(x);
    if (r < 1) return { fate: 'horizon', steps: i };
    if (r > rEsc && dot(x, v) > 0) {
      // rotate v towards the centre by the remaining bend
      const bend = tailBend(x, v);
      const vn = add([0, 0, 0], v, 1 / len(v));
      const n = cross(cross(x, vn), vn); // points roughly towards the centre, perpendicular to v
      const nl = len(n);
      const u = nl > 0 ? add([0, 0, 0], n, 1 / nl) : [0, 0, 0];
      const d = [vn[0] * Math.cos(bend) + u[0] * Math.sin(bend), vn[1] * Math.cos(bend) + u[1] * Math.sin(bend), vn[2] * Math.cos(bend) + u[2] * Math.sin(bend)];
      return { fate: 'sky', dir: d, steps: i };
    }
    const dt = k * r * (r > 8 ? Math.min(r / 8, 4) : 1) * (r < 2.2 ? 1 : 1);
    // rk4
    const a1 = accel(x, h2);
    const x2 = add(x, v, dt / 2), v2 = add(v, a1, dt / 2), a2 = accel(x2, h2);
    const x3 = add(x, v2, dt / 2), v3 = add(v, a2, dt / 2), a3 = accel(x3, h2);
    const x4 = add(x, v3, dt), v4 = add(v, a3, dt), a4 = accel(x4, h2);
    x = [0, 1, 2].map((j) => x[j] + dt / 6 * (v[j] + 2 * v2[j] + 2 * v3[j] + v4[j]));
    v = [0, 1, 2].map((j) => v[j] + dt / 6 * (a1[j] + 2 * a2[j] + 2 * a3[j] + a4[j]));
  }
  return { fate: 'lost', steps: maxSteps };
}

// a hovering camera at radius r measures angle alpha from the outward radial
// direction; in schwarzschild coordinates the ray leaves at tan(beta) = tan(alpha) / sqrt(1 - 1/r).
export function cameraRay(r0, alpha) {
  const f = Math.sqrt(1 - 1 / r0);
  return { x: [r0, 0, 0], v: [Math.cos(alpha) * f, Math.sin(alpha), 0] };
}

// ---- matter: a probe in free fall ----
// per unit mass, proper time tau, r_s = 1 (so G M = 1/2), exactly:
//   d2x/dtau2 = -x / (2 r^3) - 3/2 L^2 x / r^5,  L = |x cross dx/dtau|
//   dt/dtau   = E / (1 - 1/r)
// x is the same flat stand-in for schwarzschild r and angles the light rays use.

// speed a hovering observer at r must give something for it to circle (fraction of c)
export const circularSpeed = (r) => Math.sqrt(0.5 / (r - 1));

// thrown from radius |x0| in the local direction d (unit), at speed w (fraction of c),
// as measured by the hovering observer there
export function throwProbe(x0, d, w) {
  const r = len(x0), rh = [x0[0] / r, x0[1] / r, x0[2] / r];
  const sf = Math.sqrt(1 - 1 / r), gam = 1 / Math.sqrt(1 - w * w);
  const dr = dot(d, rh);
  const v = [0, 1, 2].map((i) => gam * w * (dr * sf * rh[i] + (d[i] - dr * rh[i])));
  return { x: x0.slice(), v, E: gam * sf, tau: 0, t: 0 };
}

function pAcc(x, L2, out) {
  const r2 = x[0] * x[0] + x[1] * x[1] + x[2] * x[2], r = Math.sqrt(r2);
  const k = -0.5 / (r2 * r) - 1.5 * L2 / (r2 * r2 * r);
  out[0] = k * x[0]; out[1] = k * x[1]; out[2] = k * x[2];
}
// advance by dt of far-away (coordinate) time. returns false if it never got there
// (only near the horizon, where its own time almost stops).
const _a = [[0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0]], _x = [0, 0, 0], _v = [0, 0, 0];
export function stepProbe(p, dt, opts = {}) {
  const kStep = opts.k ?? 0.01, maxSub = opts.maxSub ?? 4000;
  let left = dt;
  const x = p.x, v = p.v;
  for (let s = 0; s < maxSub && left > 1e-12; s++) {
    const c = [x[1] * v[2] - x[2] * v[1], x[2] * v[0] - x[0] * v[2], x[0] * v[1] - x[1] * v[0]];
    const L2 = c[0] * c[0] + c[1] * c[1] + c[2] * c[2];
    const r = len(x), f = 1 - 1 / r;
    if (r < 1 + 1e-9) return false;
    // never step more than a fifth of the way to the horizon: seen from far away it only creeps closer
    let h = Math.min(kStep * r * Math.sqrt(r), 0.2 * (r - 1));
    const hMax = left * f / p.E;   // the proper time that fits in what is left
    if (h > hMax) h = hMax;
    // rk4 on (x, v, t)
    const A = _a;
    pAcc(x, L2, A[0]);
    for (let i = 0; i < 3; i++) { _x[i] = x[i] + v[i] * h / 2; _v[i] = v[i] + A[0][i] * h / 2; }
    const r2 = len(_x); pAcc(_x, L2, A[1]); const v2 = _v.slice();
    for (let i = 0; i < 3; i++) { _x[i] = x[i] + v2[i] * h / 2; _v[i] = v[i] + A[1][i] * h / 2; }
    const r3 = len(_x); pAcc(_x, L2, A[2]); const v3 = _v.slice();
    for (let i = 0; i < 3; i++) { _x[i] = x[i] + v3[i] * h; _v[i] = v[i] + A[2][i] * h; }
    const r4 = len(_x); pAcc(_x, L2, A[3]); const v4 = _v.slice();
    const tdot = (rr) => p.E / (1 - 1 / rr);
    const dtc = h / 6 * (tdot(r) + 2 * tdot(r2) + 2 * tdot(r3) + tdot(r4));
    for (let i = 0; i < 3; i++) {
      x[i] += h / 6 * (v[i] + 2 * v2[i] + 2 * v3[i] + v4[i]);
      v[i] += h / 6 * (A[0][i] + 2 * A[1][i] + 2 * A[2][i] + A[3][i]);
    }
    p.tau += h; p.t += dtc; left -= dtc;
  }
  return left <= 1e-9;
}

// what the hovering observer at the probe's place measures
export function probeLocal(p) {
  const x = p.x, v = p.v, r = len(x), rh = [x[0] / r, x[1] / r, x[2] / r];
  const f = 1 - 1 / r, sf = Math.sqrt(f), vr = dot(v, rh);
  const u = [0, 1, 2].map((i) => vr / sf * rh[i] + (v[i] - vr * rh[i]));
  const gam = Math.sqrt(1 + dot(u, u));
  return { r, w: u.map((c) => c / gam), gamma: gam, speed: Math.sqrt(dot(u, u)) / gam, clock: f / p.E };
}

// ---- where a point appears: the light ray from the camera that passes through it ----
// returns the local direction (unit, camera frame) of the main image of P, and whether
// that light crossed the disk on the way. guess: last answer's angle, for speed.
const R0 = { r: 0 }, RINF = { r: Infinity };
export function imageOf(cam, P, guess, diskOut = DISK_OUT, tol = 1.5e-4) {
  const r0 = len(cam), ch = [cam[0] / r0, cam[1] / r0, cam[2] / r0];
  const pc = dot(P, ch);
  let e2 = [P[0] - pc * ch[0], P[1] - pc * ch[1], P[2] - pc * ch[2]];
  const pl = len(e2);
  if (pl < 1e-9) return null;
  e2 = e2.map((c) => c / pl);
  const Phi = Math.atan2(pl, pc), rP = len(P);
  const sf0 = Math.sqrt(1 - 1 / r0);
  const rEsc = Math.max(3 * rP, 1.5 * r0, 60);
  // shoot at angle psi from the inward direction, towards P's side
  const shoot = (psi) => {
    // in the plane (ch, e2): start at (r0, 0); the radial part of the direction is scaled
    let x0 = r0, x1 = 0, v0 = -Math.cos(psi) * sf0, v1 = Math.sin(psi);
    const h2 = (x0 * v1 - x1 * v0) ** 2, c1 = ch[1], e21 = e2[1];
    let hidden = false, th = 0;
    for (let i = 0; i < 600; i++) {
      const r = Math.sqrt(x0 * x0 + x1 * x1);
      if (r < 1) return R0;
      if (r > rEsc && x0 * v0 + x1 * v1 > 0) return RINF;
      const dt = STEP_K * r * Math.min(Math.max(r / 8, 1), 4) * 0.6;
      let q = x0 * x0 + x1 * x1, k = -1.5 * h2 / (q * q * Math.sqrt(q));
      const a0 = k * x0, a1 = k * x1;
      let p0 = x0 + v0 * dt / 2, p1 = x1 + v1 * dt / 2; q = p0 * p0 + p1 * p1; k = -1.5 * h2 / (q * q * Math.sqrt(q));
      const b0 = k * p0, b1 = k * p1, w0 = v0 + a0 * dt / 2, w1 = v1 + a1 * dt / 2;
      p0 = x0 + w0 * dt / 2; p1 = x1 + w1 * dt / 2; q = p0 * p0 + p1 * p1; k = -1.5 * h2 / (q * q * Math.sqrt(q));
      const c0 = k * p0, cc1 = k * p1, y0 = v0 + b0 * dt / 2, y1 = v1 + b1 * dt / 2;
      p0 = x0 + y0 * dt; p1 = x1 + y1 * dt; q = p0 * p0 + p1 * p1; k = -1.5 * h2 / (q * q * Math.sqrt(q));
      const d0 = k * p0, d1 = k * p1, z0 = v0 + c0 * dt, z1 = v1 + cc1 * dt;
      const n0 = x0 + dt / 6 * (v0 + 2 * w0 + 2 * y0 + z0), n1 = x1 + dt / 6 * (v1 + 2 * w1 + 2 * y1 + z1);
      v0 += dt / 6 * (a0 + 2 * b0 + 2 * c0 + d0);
      v1 += dt / 6 * (a1 + 2 * b1 + 2 * cc1 + d1);
      const ya = x0 * c1 + x1 * e21, yb = n0 * c1 + n1 * e21;   // height above the disk plane
      let th2 = Math.atan2(n1, n0); if (th2 < th - Math.PI) th2 += 2 * Math.PI;
      const kk = ya * yb < 0 ? ya / (ya - yb) : 2;
      if (th2 >= Phi) {
        const kc = (Phi - th) / Math.max(th2 - th, 1e-12);
        if (kk < kc) { const rd = Math.hypot(x0 + (n0 - x0) * kk, x1 + (n1 - x1) * kk); if (rd > 3 && rd < diskOut) hidden = true; }
        return { r: Math.hypot(x0 + (n0 - x0) * kc, x1 + (n1 - x1) * kc), hidden };
      }
      if (kk <= 1) { const rd = Math.hypot(x0 + (n0 - x0) * kk, x1 + (n1 - x1) * kk); if (rd > 3 && rd < diskOut) hidden = true; }
      x0 = n0; x1 = n1; th = th2;
    }
    return R0;
  };
  const F = (psi) => { const s = shoot(psi); return { s, g: s.r === Infinity ? 1e9 : s.r - rP }; };
  let lo = 0, hi = Math.PI, best = null;
  if (guess != null) {
    const a = Math.max(1e-4, guess - 0.006), b = Math.min(Math.PI - 1e-4, guess + 0.006);
    const fa = F(a), fb = F(b);
    if (fa.g < 0 && fb.g > 0) { lo = a; hi = b; }
  }
  for (let i = 0; i < 40 && hi - lo > tol; i++) {
    const m = (lo + hi) / 2, fm = F(m);
    if (fm.g < 0) lo = m; else { hi = m; best = fm.s; }
  }
  const psi = (lo + hi) / 2;
  if (!best || best.r === Infinity) return null;
  const d = [0, 1, 2].map((i) => -Math.cos(psi) * ch[i] + Math.sin(psi) * e2[i]);
  return { d, psi, hidden: best.hidden };
}
