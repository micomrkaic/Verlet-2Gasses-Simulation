# Parabolic Vessel — Hard Sphere Gas Simulation

A real-time 2D molecular dynamics simulation of hard-sphere particles in a
configurable vessel, written in C with SDL2. Demonstrates thermalisation,
energy conservation, Maxwell-Boltzmann statistics, equipartition between
species, and acoustic compression waves via an oscillating piston ceiling.
Supports video export to MP4.

**[▶ Play in browser](https://micomrkaic.itch.io/two-d-gas-box-with-verlet)**
(WebAssembly build on itch.io)

---

## Two builds

| | v1 | v2 |
|---|---|---|
| Source | `verlet_collisions_two_gasses_video.c` | `ballsim_v2.c` |
| Binary | `ballsim` | `ballsim_v2` |
| Algorithm | Event-driven MD, exact CCD | Fixed-dt Verlet + cell-list |
| Max particles | 500 | 5 000 |
| Collision detection | O(N²) pairwise | O(N) via cell list |
| Parallelism | None | OpenMP (Verlet step) |
| Energy conservation | < 0.15 % drift | < 1 % typical |
| Best for | Exact physics, small N | High particle counts |

Both binaries share identical controls, sliders, stats panel, and video
recording — only the physics engine differs.

---

## Physics

Particles are hard spheres under constant downward gravity inside a vessel
whose floor is a parabolic curve:

```
y_floor(x) = PARA_BASE - para_a × (x − sim_w/2)²
```

Setting `para_a = 0` gives a flat rectangular box. Vertical side walls and a
movable ceiling piston complete the enclosure.

**Collisions are fully elastic** by default (`restitution = 1.0`).

### v1 — event-driven MD (exact)

- **Particle–parabola / walls**: contact time solved analytically as a
  quadratic equation. The trajectory under gravity is substituted into the
  parabola equation, yielding exact coefficients `A`, `B`, `C` whose roots
  give the precise contact moment. Velocity is reflected about the surface
  normal.
- **Particle–particle**: fully event-driven MD — all pairwise contact times
  are solved analytically (gravity cancels in the relative frame, making the
  quadratic exact). The earliest event wins; all particles advance to that
  moment and exchange velocities via the elastic impulse formula.
- **Ceiling piston**: binary search finds the exact contact time with the
  moving surface; the bounce reflects velocity relative to the piston.

Energy conservation is **< 0.15% drift over 12 000 frames** with a mixed
two-species gas.

### v2 — cell-list Verlet MD

Each outer timestep (1/120 s) is divided into **8 sub-steps**:

1. **Verlet integrate** under gravity — O(N), OpenMP-parallel
2. **Boundary conditions** — walls, ceiling, parabola floor pushed out and
   reflected analytically
3. **Build cell list** — O(N); cell size = 2 × max_radius so only the
   3 × 3 = 9 neighbouring cells need checking per particle
4. **Resolve overlaps** — 4 passes of impulse + position correction per
   sub-step; handles dense packings without drift accumulation

The cell list reduces collision detection from O(N²) to O(N), enabling
thousands of particles at 60 fps. OpenMP parallelises the Verlet step across
all available cores.

### Two-species gas

Particles come in two species — **light** (cyan) and **heavy** (orange):

| Property | Light | Heavy |
|----------|-------|-------|
| Mass | 1 | `mass_ratio` |
| Radius | `radius` | `radius × ∛mass_ratio` |
| Color | cyan-blue | orange-red |

Radius scales as the cube root of mass (equal density, volume ∝ mass).

Elastic collisions between unequal masses use the correct impulse formula:

```
J = 2·mᵢ·mⱼ / (mᵢ + mⱼ) · dvₙ
```

where `dvₙ` is the relative normal velocity at contact.

### Equipartition of energy

At equilibrium the equipartition theorem requires each species to have the
same mean KE per particle:

```
½·m_light·⟨v_light²⟩  =  ½·m_heavy·⟨v_heavy²⟩  =  kT
```

The LIGHT KE/N and HEAVY KE/N panel readouts let you watch convergence in
real time.

### Thermalisation

Starting monoenergetic (all same KE), collisions redistribute energy and the
distribution converges to the 2D Maxwell-Boltzmann exponential:

```
f(KE) ∝ exp(−KE / kT)    where kT = mean KE per particle
```

### Oscillating piston

When `ceil_amp > 0` the ceiling acts as a sinusoidal piston:

```
y_ceil(t) = radius + ceil_amp × (1 + sin(2π × ceil_freq × t))
```

For near-zero net energy drift the piston must be quasi-static:

```
ceil_amp × 2π × ceil_freq  ≪  mono_speed   (ratio < ~0.1)
```

When faster, irreversible **Fermi acceleration** heats the gas — real
physics, not a bug.

---

## Building

### Dependencies

- **SDL2** — display and input
- **ffmpeg** — video export (optional)
- **GCC with OpenMP** — for v2 parallel Verlet (falls back to
  single-threaded without `-fopenmp`)

```bash
# Linux
apt install libsdl2-dev ffmpeg

# macOS (Homebrew)
brew install sdl2 ffmpeg
```

### Quick build (both binaries)

```bash
make          # builds ./ballsim (v1) and ./ballsim_v2 (v2)
make run      # launch v1
make run2     # launch v2
make clean
```

### Manual build

```bash
# v1 — event-driven
gcc -O2 -std=c11 -Wall $(pkg-config --cflags sdl2) \
    -o ballsim verlet_collisions_two_gasses_video.c \
    $(pkg-config --libs sdl2) -lm

# v2 — cell-list + OpenMP
gcc -O3 -std=c11 -Wall -fopenmp $(pkg-config --cflags sdl2) \
    -o ballsim_v2 ballsim_v2.c \
    $(pkg-config --libs sdl2) -lm
```

### WebAssembly (itch.io / browser)

Requires [Emscripten](https://emscripten.org/docs/getting_started/downloads.html):

```bash
source /path/to/emsdk/emsdk_env.sh
make -f Makefile.emscripten          # → ballsim.html + .js + .wasm
python3 -m http.server 8080          # test locally
# open http://localhost:8080/ballsim.html
```

To deploy on itch.io:

```bash
zip ballsim.zip ballsim.html ballsim.js ballsim.wasm
```

Upload the ZIP to itch.io, tick **"This file will be played in the
browser"**, set embed dimensions to **980 × 960**, and tick **"Enable
SharedArrayBuffer support"**.

---

## Configuration

Parameters load from a plain-text config file at startup. Most can also be
adjusted live via the **in-app sliders**.

```bash
./ballsim    --write-config ballsim.cfg    # write default template
./ballsim_v2 --write-config ballsim_v2.cfg

$EDITOR ballsim.cfg
./ballsim --config ballsim.cfg
```

### Box dimensions

`sim_w` and `sim_h` require a window resize — set them in the config file
and restart. They cannot be changed via sliders.

### Parameter reference

| Parameter | v1 default | v2 default | Slider | Effect |
|-----------|-----------|-----------|--------|--------|
| `sim_w` | 700 | 700 | — | Simulation area width (restart required) |
| `sim_h` | 720 | 960 | — | Simulation area height (restart required) |
| `n_particles` | 150 | 300 | ✓ | Starting particle count |
| `gravity` | 500 px/s² | 500 px/s² | ✓ | Gravitational acceleration |
| `radius` | 5 px | 4 px | ✓ | Base radius for light particles |
| `restitution` | 1.0 | 1.0 | ✓ | < 1.0 adds damping |
| `mono_speed` | 220 px/s | 220 px/s | ✓ | Speed for M-key reset |
| `para_a` | 0.003 | 0.003 | ✓ | Parabola curvature; 0 = flat box |
| `ceil_amp` | 0 px | 0 px | ✓ | Piston amplitude; 0 = fixed ceiling |
| `ceil_freq` | 0.5 Hz | 0.5 Hz | ✓ | Piston oscillation frequency |
| `heavy_frac` | 0.5 | 0.5 | ✓ | Fraction of heavy particles (0–1) |
| `mass_ratio` | 4.0 | 4.0 | ✓ | Mass of heavy / mass of light (1–20) |

**Particle slider range:** 1–500 in v1, 1–5000 in v2.

**v2 guidance for high N:** reduce `radius` as you increase `n_particles`
to keep packing fraction reasonable. At N = 2000 with `radius = 3` and
`mono_speed = 150` the simulation runs comfortably at 60 fps on a 4-core
machine. At N = 5000 use `radius = 2`.

---

## In-app sliders

The lower portion of the right-hand panel contains ten sliders. Click
anywhere on a track to jump to that value; hold and drag to scrub
continuously. The active slider highlights in bright cyan. Changes that
affect particle placement (radius, para_a, n_particles, mono_speed,
heavy_frac, mass_ratio) trigger an automatic reset.

---

## Controls

| Key | Action |
|-----|--------|
| `M` | Reset — monoenergetic (equal KE per particle) |
| `R` | Reset — random speeds |
| `V` | Toggle video recording (v2 panel shows OMP thread count instead of V-REC) |
| `+` / `=` | Add one particle |
| `-` | Remove one particle |
| `Q` / `Esc` | Quit |

---

## Stats panel

| Readout | Meaning |
|---------|---------|
| BOX SIZE | sim_w × sim_h in pixels |
| PARTICLES | Live particle count |
| TOTAL ENERGY (E/E0) | Energy normalised to initial — 1.0000 = conserved |
| KINETIC (KE/E) | Kinetic energy as fraction of total |
| POTENTIAL (PE/E) | Potential energy as fraction of total |
| TEMPERATURE (KE/N) | Mean KE per particle across all species |
| LIGHT KE/N | Mean KE per light particle |
| HEAVY KE/N | Mean KE per heavy particle |
| COLLISIONS | Cumulative particle–particle collision count |
| V2.0 badge | Shown in v2 only |
| OMP THREADS-N | OpenMP thread count (v2 only) |

### Energy trace

Ring buffer of the last ~280 samples of E/E0. White line = 1.0 (perfect
conservation). Bars colour from green toward red as drift grows. DRIFT
shows the current percentage deviation.

### KE distribution

Histogram of per-particle KE (cyan bars, exponentially smoothed) overlaid
with the theoretical 2D Maxwell-Boltzmann curve (yellow). Press `M` to reset
to a monoenergetic spike and watch it relax to the exponential shape.

---

## Video recording

Press `V` to start; `V` or `Q` to stop. Output:

```
ballsim_YYYYMMDD_HHMMSS.mp4     (v1)
ballsim_v2_YYYYMMDD_HHMMSS.mp4  (v2)
```

**Format:** H.264, YUV420p, 60 fps, near-lossless (`-crf 18`). Frames are
piped directly into ffmpeg — no temporary files. Capture is gated on
wall-clock time so playback speed always matches what you see on screen.

Requires `ffmpeg` on your `PATH`. Not available in the WASM/itch.io build.

---

## Repository layout

```
.
├── verlet_collisions_two_gasses_video.c   ← v1 source (event-driven, untouched)
├── ballsim_v2.c                           ← v2 source (cell-list MD + OpenMP)
├── ballsim_wasm.c                         ← WASM source (v1 physics, browser build)
├── Makefile                               ← builds v1 + v2 natively
├── Makefile.emscripten                    ← builds WASM with emcc
├── shell.html                             ← custom HTML wrapper for WASM build
└── README.md
```

---

## Suggested experiments

**Equipartition** — press `M` with `mass_ratio = 8`, `heavy_frac = 0.5`.
Watch LIGHT KE/N and HEAVY KE/N converge from equal initial values to the
same equilibrium temperature.

**Mass segregation** — high `mass_ratio` (try 16) + high `gravity`: heavy
particles settle near the bowl bottom, light ones bounce higher.

**Thermalisation** — press `M`, watch the KE histogram spike broaden into
the exponential MB curve over a few seconds.

**Flat box** — set `para_a = 0` for a rectangular box with uniform gravity.

**Compression waves** — `para_a = 0`, `ceil_amp = 5`, `ceil_freq = 0.5`.
Density waves propagate downward on each piston downstroke.

**Inelastic collapse** — set `restitution = 0.8`, watch kinetic energy drain
slowly from the energy trace.

**Dense gas (v2)** — in v2, set N = 2000, radius = 3. The collision rate is
high enough that thermalisation is nearly instantaneous and the MB curve is
always sharp.

**Scale to thousands (v2 only)** — drag the PARTICLES slider toward 5000
while reducing radius. The cell-list engine keeps 60 fps up to the packing
limit. OpenMP thread count is shown at the bottom of the panel.
