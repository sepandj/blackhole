// node test/physics.test.mjs
// checks the ray integrator and the disk model against textbook results.
import * as P from '../physics.js';

let fails = 0;
function check(name, got, want, tol) {
  const err = Math.abs(got - want) / Math.abs(want);
  const ok = err <= tol;
  if (!ok) fails++;
  console.log(`${ok ? 'ok  ' : 'FAIL'} ${name}: got ${got.toPrecision(7)} want ${want.toPrecision(7)} (rel err ${err.toExponential(1)})`);
}

// a ray from far away with impact parameter b, aimed past the hole
function farRay(b, r0 = 2000) {
  return { x: [r0, 0, 0], v: [-Math.sqrt(1 - (b / r0) ** 2), b / r0, 0] };
}
function deflection(b) {
  const { x, v } = farRay(b);
  const t = P.trace(x, v, { maxSteps: 4000 });
  if (t.fate !== 'sky') return NaN;
  // started heading in -x; total bend is the angle between start and end directions
  const d = t.dir; const v0 = v;
  // signed angle in the orbit plane, unwrapped (rays near the photon sphere bend past pi)
  const a0 = Math.atan2(v0[1], v0[0]), a1 = Math.atan2(d[1], d[0]);
  let da = a1 - a0; while (da < 0) da += 2 * Math.PI;
  return da;
}

// exact deflection for schwarzschild (r_s = 1): from (du/dphi)^2 = 1/b^2 - u^2 + u^3
function exactDeflection(b) {
  // turning point u0: largest root below 2/3 of 1/b^2 - u^2 + u^3 = 0
  const g = (u) => 1 / (b * b) - u * u + u * u * u;
  let lo = 0, hi = 2 / 3;
  for (let i = 0; i < 200; i++) { const m = (lo + hi) / 2; if (g(m) > 0) lo = m; else hi = m; }
  const u0 = lo;
  // substitute u = u0 (1 - s^2) to remove the endpoint singularity
  const n = 200000; let sum = 0;
  for (let i = 0; i < n; i++) {
    const s = (i + 0.5) / n;
    const u = u0 * (1 - s * s);
    sum += (2 * u0 * s) / Math.sqrt(Math.max(g(u), 1e-300));
  }
  return 2 * sum / n - Math.PI;
}

console.log('-- light bending --');
for (const b of [100, 30]) {
  const weak = 2 / b + 15 * Math.PI / (16 * b * b);
  check(`deflection b=${b} r_s vs weak-field series`, deflection(b), weak, 2e-3 * (100 / b));
}
for (const b of [10, 6, 4, 3, 2.75, 2.65]) {
  check(`deflection b=${b} r_s vs exact integral`, deflection(b), exactDeflection(b), 2e-3);
}
// the sun, with its real mass and radius: 1.75 arcseconds
{
  const rs = P.schwarzschildRadius(P.MSUN), b = 6.957e8 / rs;
  const arcsec = (2 / b) * 180 / Math.PI * 3600;
  check('light grazing the sun (arcsec)', arcsec, 1.751, 2e-3);
}

console.log('-- shadow --');
{
  let lo = 2, hi = 3.5;
  for (let i = 0; i < 40; i++) {
    const m = (lo + hi) / 2;
    const { x, v } = farRay(m);
    const t = P.trace(x, v, { maxSteps: 6000 });
    if (t.fate === 'horizon') lo = m; else hi = m;
  }
  check('critical impact parameter (sqrt 27 / 2)', (lo + hi) / 2, P.SHADOW_B, 1e-3);
}
// a hovering camera close in: synge's formula sin(a) = (sqrt 27 / 2) / r * sqrt(1 - 1/r)
for (const r0 of [20, 6, 3, 1.8]) {
  let lo = 0, hi = Math.PI; // angle from the outward direction; inward is pi
  // shadow is around the inward direction: rays with alpha near pi fall in
  for (let i = 0; i < 40; i++) {
    const m = (lo + hi) / 2;
    const { x, v } = P.cameraRay(r0, m);
    const t = P.trace(x, v, { maxSteps: 6000 });
    if (t.fate === 'horizon') hi = m; else lo = m;
  }
  const shadowHalf = Math.PI - (lo + hi) / 2;
  const want = r0 > 1.5 ? Math.asin(P.SHADOW_B / r0 * Math.sqrt(1 - 1 / r0)) : Math.PI - Math.asin(P.SHADOW_B / r0 * Math.sqrt(1 - 1 / r0));
  check(`shadow half-angle from r=${r0} r_s (synge)`, shadowHalf, want, 2e-3);
}

console.log('-- disk --');
check('efficiency 1 - sqrt(8/9)', P.ETA, 0.0571909584, 1e-9);
{
  // page & thorne 1974 closed form for a = 0, x = sqrt(r / M), x0 = sqrt 6
  const s3 = Math.sqrt(3), x0 = Math.sqrt(6);
  const closed = (r) => {
    const x = Math.sqrt(r);
    const br = x - x0 - s3 / 2 * Math.log((x - s3) / (x0 - s3)) + s3 / 2 * Math.log((x + s3) / (x0 + s3));
    return 3 / (8 * Math.PI) * br / (x ** 4 * (x ** 3 - 3 * x));
  };
  for (const r of [6.5, 9.55, 20, 100, 2000]) check(`page-thorne flux at r=${r} M vs closed form`, P.pageThorneFluxM(r), closed(r), 2e-4);
  const pk = P.diskPeak();
  check('flux peak radius (r_s), known ~4.77', pk.r, 4.77, 1e-2);
}
{
  // luminosity of the disk: both faces integrated with redshift should be eta * mdot c^2.
  // at infinity: L = 2 * integral F * E * 2 pi r dr  (page-thorne, geometric units)
  let L = 0; const dr = 0.01;
  for (let r = 6 + dr / 2; r < 40000; r += dr * Math.max(1, r / 50)) {
    const s = Math.sqrt(1 - 3 / r), E = (1 - 2 / r) / s;
    const step = dr * Math.max(1, r / 50);
    L += 2 * P.pageThorneFluxM(r) * E * 2 * Math.PI * r * step;
  }
  check('disk luminosity / mdot c^2 = eta', L, P.ETA, 2e-2);
}

console.log('-- colour --');
{
  const sun = P.blackbody(5772);
  console.log('sun rgb', sun.rgb.map((x) => x.toFixed(3)).join(' '));
  check('sun-like blackbody surface luminance (cd/m^2), ~2.0e9', sun.lum, 2.0e9, 0.1);
}

console.log('-- the three holes --');
for (const h of P.HOLES) {
  const m = h.mass * P.MSUN, rs = P.schwarzschildRadius(m);
  const pk = P.diskPeak();
  for (const e of [1e-1, 1e-3, 1e-5]) {
    const T = P.diskT0(m, e) * Math.pow(pk.f, 0.25);
    console.log(`${h.name.padEnd(16)} r_s ${(rs / 1e3).toExponential(3)} km  edd ${e}  peak T ${T.toExponential(2)} K  isco period ${(2 * Math.PI / P.orbitOmega(3) * P.timeUnit(m) / 3600).toPrecision(3)} h`);
  }
}

console.log(fails ? `\n${fails} failed` : '\nall passed');
process.exit(fails ? 1 : 0);
