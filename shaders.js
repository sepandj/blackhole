// glsl for the black hole. units: schwarzschild radius r_s = 1, time in r_s / c.
// the disk lies in the x-z plane and turns about +y.

export const VERT = `#version 300 es
in vec2 aPos;
out vec2 vUv;
void main(){ vUv = aPos * 0.5 + 0.5; gl_Position = vec4(aPos, 0.0, 1.0); }`;

export const TRACE = `#version 300 es
precision highp float;
precision highp sampler2D;
uniform vec2 uRes;
uniform vec3 uCam;        // camera position
uniform mat3 uBasis;      // right, up, forward
uniform vec2 uTan;        // tan of half field of view, x and y
uniform float uF;         // sqrt(1 - 1/r0): the camera's clock rate
uniform float uTime;      // coordinate time
uniform float uT0;        // disk temperature scale, kelvin
uniform float uExpo;      // display units per cd/m^2
uniform float uSky;       // sky brightness
uniform vec3 uSkyGain;    // blueshift of a 5000 K sky seen by the hovering camera
uniform float uPix;       // angle of one pixel, radians
uniform int uSteps;
uniform float uDiskOut;
uniform sampler2D uDiskT; // F^(1/4) over the disk, log radius
uniform sampler2D uBB;    // blackbody: rgb chromaticity, log10 luminance
out vec4 o;

const float K = 0.075;    // step = K r, same as physics.js
const float BB_LO = 300.0, BB_HI = 1e9;

vec4 bbRaw(float T){
  float u = clamp(log(T / BB_LO) / log(BB_HI / BB_LO), 0.0, 1.0) * 511.0;
  int i = min(int(u), 510);
  return mix(texelFetch(uBB, ivec2(i, 0), 0), texelFetch(uBB, ivec2(i + 1, 0), 0), u - float(i));
}
vec3 bb(float T){ vec4 t = bbRaw(T); return t.rgb * pow(10.0, t.a); }

float hash13(vec3 p){ p = fract(p * 0.1031); p += dot(p, p.zyx + 31.32); return fract((p.x + p.y) * p.z); }
vec3 hash33(vec3 p){
  p = fract(p * vec3(0.1031, 0.1030, 0.0973));
  p += dot(p, p.yxz + 33.33);
  return fract((p.xxy + p.yxx) * p.zyx);
}
float vnoise(vec3 p){
  vec3 i = floor(p), f = fract(p);
  vec3 u = f * f * (3.0 - 2.0 * f);
  float a = hash13(i), b = hash13(i + vec3(1,0,0)), c = hash13(i + vec3(0,1,0)), d = hash13(i + vec3(1,1,0));
  float e = hash13(i + vec3(0,0,1)), f1 = hash13(i + vec3(1,0,1)), g = hash13(i + vec3(0,1,1)), h = hash13(i + vec3(1,1,1));
  return mix(mix(mix(a,b,u.x), mix(c,d,u.x), u.y), mix(mix(e,f1,u.x), mix(g,h,u.x), u.y), u.z);
}
float fbm(vec3 p){
  float s = 0.0, a = 0.5;
  for (int i = 0; i < 5; i++){ s += a * vnoise(p); p = p * 2.03 + 17.1; a *= 0.5; }
  return s;
}

// ---- the sky: stars with real colours, and a galaxy band ----
const vec3 GAL_N = vec3(0.276, 0.906, -0.320);   // normal of the galactic plane (tilted from the disk)
const vec3 GAL_C = vec3(0.857, -0.188, 0.480);   // towards the galactic centre (in the plane)

vec3 starLayer(vec3 s, float n, float density, float bright, float seed){
  vec3 q = s * n;
  vec3 c = floor(q);
  vec3 h = hash33(c + seed);
  if (h.x > density) return vec3(0);
  vec3 sp = normalize(c + 0.2 + 0.6 * hash33(c + seed + 7.0));
  float d = length(s - sp);
  float sig = uPix * 0.75;
  float flux = bright * pow(h.y, 14.0) + bright * 0.006 * h.y;
  // temperature: most stars are cool and orange, a few hot and blue
  float T = 2800.0 + 9000.0 * pow(h.z, 2.2) + 22000.0 * pow(h.z, 30.0);
  vec4 em = bbRaw(T), ob = bbRaw(T / uF);
  float gain = pow(10.0, ob.a - em.a);
  return ob.rgb * flux * gain * exp(-0.5 * d * d / (sig * sig)) / (sig * sig) * 1e-6;
}

vec3 sky(vec3 s){
  vec3 col = vec3(0);
  col += starLayer(s, 34.0, 0.22, 22.0, 1.0);
  col += starLayer(s, 80.0, 0.14, 5.0, 2.0);
  col += starLayer(s, 190.0, 0.10, 1.6, 3.0);
  // galaxy band, brighter towards its centre, cut by dust lanes
  float lat = dot(s, GAL_N);
  float toC = dot(s, GAL_C) * 0.5 + 0.5;
  float band = exp(-lat * lat / 0.018) * (0.25 + 0.75 * pow(toC, 3.0));
  float bulge = exp(-(lat * lat) / 0.06) * pow(toC, 14.0);
  float cloud = fbm(s * 4.0 + 3.0);
  float dust = smoothstep(0.42, 0.75, fbm(s * 9.0 + vec3(5.0, 1.0, 2.0) + lat * 6.0));
  float dustLane = exp(-pow((lat - 0.012) / 0.03, 2.0)) * smoothstep(0.35, 0.7, fbm(s * 14.0));
  float glow = (band * (0.45 + 0.9 * cloud) + bulge * 1.4) * (1.0 - 0.75 * dust) * (1.0 - 0.8 * dustLane);
  vec3 tint = mix(vec3(0.55, 0.62, 0.85), vec3(1.0, 0.82, 0.62), clamp(toC * toC + bulge, 0.0, 1.0));
  col += tint * glow * 0.035;
  // a faint grainy unresolved star field so the band has texture
  col += tint * band * 0.05 * pow(hash13(floor(s * 900.0)), 18.0);
  return col * uSkyGain * uSky;
}

// ---- the disk ----
float diskTau(float r){
  float u = clamp(log(r / 3.0) / log(uDiskOut / 3.0), 0.0, 1.0) * 255.0;
  int i = min(int(u), 254);
  return mix(texelFetch(uDiskT, ivec2(i, 0), 0).r, texelFetch(uDiskT, ivec2(i + 1, 0), 0).r, u - float(i));
}

// turbulence that orbits with the gas: each ring turns at its own kepler rate.
// two layers fade into each other so the shear never winds up forever.
float diskDensity(float r, float phi, float om){
  const float P = 260.0;
  float t1 = mod(uTime, P), t2 = mod(uTime + 0.5 * P, P);
  float w = abs(1.0 - 2.0 * t1 / P);
  float k1 = floor(uTime / P), k2 = floor((uTime + 0.5 * P) / P);
  float lr = log(r);
  float a1 = phi + om * t1, a2 = phi + om * t2;
  // sampled in real distance so detail is even across the disk; stretched along the orbit
  float n1 = fbm(vec3(cos(a1) * r * 0.55, sin(a1) * r * 0.55, r * 2.4) + k1 * 13.7);
  float n2 = fbm(vec3(cos(a2) * r * 0.55, sin(a2) * r * 0.55, r * 2.4) + k2 * 13.7);
  float n = mix(n2, n1, w);
  float rings = 0.5 + 0.5 * sin(lr * 90.0 + 7.0 * n);
  float streak = smoothstep(0.28, 0.72, n);
  return streak * (0.8 + 0.4 * rings);
}

vec3 acc(vec3 x, float h2){ float r2 = dot(x, x); return -1.5 * h2 * x / (r2 * r2 * sqrt(r2)); }

void main(){
  vec2 p = (gl_FragCoord.xy - 0.5 * uRes) / (0.5 * uRes);
  vec3 d = normalize(uBasis * vec3(p.x * uTan.x, p.y * uTan.y, 1.0));

  vec3 x = uCam;
  float r0 = length(x);
  vec3 rh = x / r0;
  // a hovering observer's angles map to coordinate directions with the radial part scaled by sqrt(1 - 1/r)
  vec3 v = d + dot(d, rh) * (uF - 1.0) * rh;
  vec3 Lv = cross(x, v);
  float h2 = dot(Lv, Lv);
  // conserved impact parameter: 1/b^2 = |v|^2/h^2 - 1/r^3
  float b = inversesqrt(max(dot(v, v) / max(h2, 1e-20) - 1.0 / (r0 * r0 * r0), 1e-20));
  float lamT = b * Lv.y * inversesqrt(max(h2, 1e-20));   // traced photon's L_y / E

  vec3 col = vec3(0);
  float tr = 1.0;
  float rEsc = max(1.1 * r0, 40.0);

  for (int i = 0; i < 1000; i++){
    if (i >= uSteps) { tr = 0.0; break; }
    float r = length(x);
    if (r < 1.0) { tr = 0.0; break; }
    if (r > rEsc && dot(x, v) > 0.0){
      // weak field from here on: add the remaining bend analytically
      float vl = length(v);
      vec3 vn = v / vl;
      float cps = dot(x, vn) / r;
      float bb2 = length(cross(x, vn));
      float bend = (1.0 - cps * (3.0 - cps * cps) * 0.5) / max(bb2, 1e-6);
      vec3 n = cross(cross(x, vn), vn);
      float nl = length(n);
      vec3 s = nl > 1e-9 ? vn * cos(bend) + (n / nl) * sin(bend) : vn;
      col += tr * sky(s);
      tr = 0.0;
      break;
    }
    float dt = K * r * clamp(r / 8.0, 1.0, 4.0);
    vec3 a1 = acc(x, h2);
    vec3 x2 = x + v * (0.5 * dt), v2 = v + a1 * (0.5 * dt); vec3 a2 = acc(x2, h2);
    vec3 x3 = x + v2 * (0.5 * dt), v3 = v + a2 * (0.5 * dt); vec3 a3 = acc(x3, h2);
    vec3 x4 = x + v3 * dt, v4 = v + a3 * dt; vec3 a4 = acc(x4, h2);
    vec3 xn = x + (dt / 6.0) * (v + 2.0 * v2 + 2.0 * v3 + v4);
    vec3 vn = v + (dt / 6.0) * (a1 + 2.0 * a2 + 2.0 * a3 + a4);

    if (x.y * xn.y < 0.0){
      float t = x.y / (x.y - xn.y);
      vec3 hp = mix(x, xn, t);
      float rr = length(hp);
      if (rr > 3.0 && rr < uDiskOut){
        float om = sqrt(0.5 / (rr * rr * rr));
        // redshift of gas on a circular orbit seen at infinity, then blueshift down to the hovering camera
        float g = sqrt(1.0 - 1.5 / rr) / (1.0 + om * lamT) / uF;
        float T = g * uT0 * diskTau(rr);
        float phi = atan(hp.z, hp.x);
        float dens = diskDensity(rr, phi, om);
        float edge = smoothstep(3.0, 3.25, rr) * (1.0 - smoothstep(uDiskOut * 0.55, uDiskOut, rr));
        // thick and opaque close in, thinner further out; gaps between streaks let light through
        float thick = mix(1.0, 0.3, smoothstep(5.0, 24.0, rr));
        float alpha = clamp(mix(0.12, 1.0, dens) * thick, 0.0, 1.0) * edge;
        col += tr * alpha * bb(T) * (0.45 + 0.9 * dens) * uExpo;
        tr *= 1.0 - alpha;
        if (tr < 0.01) { tr = 0.0; break; }
      }
    }
    x = xn; v = vn;
  }
  o = vec4(min(col, vec3(60000.0)), 1.0);
}`;

// dual-filter (kawase) bloom
export const DOWN = `#version 300 es
precision highp float;
uniform sampler2D uSrc; uniform vec2 uTexel;
in vec2 vUv; out vec4 o;
void main(){
  vec2 h = uTexel * 0.5;
  vec3 s = texture(uSrc, vUv).rgb * 4.0;
  s += texture(uSrc, vUv - h).rgb + texture(uSrc, vUv + h).rgb;
  s += texture(uSrc, vUv + vec2(h.x, -h.y)).rgb + texture(uSrc, vUv - vec2(h.x, -h.y)).rgb;
  o = vec4(min(s / 8.0, vec3(60000.0)), 1.0);
}`;

export const UP = `#version 300 es
precision highp float;
uniform sampler2D uSrc; uniform vec2 uTexel;
in vec2 vUv; out vec4 o;
void main(){
  vec2 h = uTexel * 0.5;
  vec3 s = texture(uSrc, vUv + vec2(-h.x * 2.0, 0.0)).rgb;
  s += texture(uSrc, vUv + vec2(-h.x, h.y)).rgb * 2.0;
  s += texture(uSrc, vUv + vec2(0.0, h.y * 2.0)).rgb;
  s += texture(uSrc, vUv + vec2(h.x, h.y)).rgb * 2.0;
  s += texture(uSrc, vUv + vec2(h.x * 2.0, 0.0)).rgb;
  s += texture(uSrc, vUv + vec2(h.x, -h.y)).rgb * 2.0;
  s += texture(uSrc, vUv + vec2(0.0, -h.y * 2.0)).rgb;
  s += texture(uSrc, vUv + vec2(-h.x, -h.y)).rgb * 2.0;
  o = vec4(s / 12.0, 1.0);
}`;

export const FINAL = `#version 300 es
precision highp float;
uniform sampler2D uScene, uBloom;
uniform float uBloomK, uFade;
uniform vec2 uRes;
in vec2 vUv; out vec4 o;
vec3 aces(vec3 x){ return clamp((x * (2.51 * x + 0.03)) / (x * (2.43 * x + 0.59) + 0.14), 0.0, 1.0); }
// compress brightness but keep each colour's hue; only the very brightest light goes to white
vec3 tone(vec3 c){
  float y = dot(c, vec3(0.2126, 0.7152, 0.0722));
  if (y <= 0.0) return vec3(0);
  float yt = aces(vec3(y)).x;
  vec3 hue = c * (yt / y);
  float m = max(hue.r, max(hue.g, hue.b));
  if (m > 1.0) hue /= m;
  return mix(hue, aces(c), smoothstep(0.6, 1.0, yt) * 0.6);
}
float ign(vec2 p){ return fract(52.9829189 * fract(dot(p, vec2(0.06711056, 0.00583715)))); }
vec3 srgb(vec3 c){ return mix(c * 12.92, 1.055 * pow(c, vec3(1.0 / 2.4)) - 0.055, step(0.0031308, c)); }
void main(){
  vec3 c = texture(uScene, vUv).rgb + uBloomK * texture(uBloom, vUv).rgb;
  c = srgb(tone(c)) * uFade;
  c += (ign(gl_FragCoord.xy) - 0.5) / 255.0;
  o = vec4(c, 1.0);
}`;
