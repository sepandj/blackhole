# black hole

https://blackhole.sepand.dev

a black hole worked out, not drawn. every pixel traces a ray of light backwards
through schwarzschild spacetime on your gpu (webgl2), using the exact photon
equation d²x/dλ² = -3/2 h² x / r⁵ (units of the horizon radius). the disk is a
novikov-thorne / page-thorne thin disk with real temperatures for real masses,
seen through its doppler and gravitational shift, coloured by the blackbody
spectrum as the eye would see it.

drag to fly around it. pinch or scroll to go closer.

- `physics.js` constants, masses, disk model, blackbody colour, the reference integrator
- `shaders.js` the ray tracer, bloom and tone mapping
- `render.js` camera, input, adaptive resolution, readouts
- `test/physics.test.mjs` checks against textbook results: `node test/physics.test.mjs`

the tests compare light bending with the exact integral (to ~10⁻⁵), the shadow
with synge's formula from far away down to 1.8 horizon radii, the disk flux with
page and thorne's closed form, its efficiency (5.7%), and 1.75 arcseconds for
light grazing the sun.

`?still` draws one full resolution frame and stops. `d`, `p`, `y` set distance
(horizon radii), pitch and yaw (degrees), `h` the hole (m87, sgra, cyg), `e` the
feeding rate as log10 of the eddington limit.
