# Parabolic Vessel — Hard Sphere Gas Simulation

A real-time 2D molecular dynamics simulation of hard-sphere particles in a
configurable vessel, written in C with SDL2. Demonstrates thermalisation,
energy conservation, Maxwell-Boltzmann statistics, equipartition between
species, barometric density profiles, the Boltzmann H-theorem, hard-disc
equations of state, and piston work via an oscillating ceiling. Supports
video export to MP4.

**[▶ Play in browser](https://micomrkaic.itch.io/two-d-gas-box-with-verlet)**
(WebAssembly build of v1 on itch.io)

---

## Two builds

| | v1 | v2.1 |
|---|---|---|
| Source | `verlet_collisions_two_gasses_video.c` | `ballsim_v2.c` |
| Binary | `ballsim` | `ballsim_v2` |
| Algorithm | Event-driven MD, exact CCD | Fixed-dt Verlet + cell list |
| Max particles | 500 | 5 000 |
| Collision detection | O(N²) pairwise | O(N) via cell list |
| Parallelism | None | OpenMP, 9-colour checkerboard |
| Rendering | Per-pixel | One batched sprite draw call |
| Energy drift (measured) | < 0.15 % / 12 000 frames | +0.002 % / 40 s at N = 1000 |
| Analysis panel | — | Yes (see below) |
| Best for | Exact physics, small N | High N, quantitative demos |

Both share the same controls, sliders, stats panel and video recording.

---

## Physics

Particles are hard spheres under constant downward gravity inside a vessel
whose floor is a parabolic curve:

```
y_floor(x) = PARA_BASE - para_a × (x − sim_w/2)²
```

Setting `para_a = 0` gives a flat rectangular box. Vertical side walls and a
movable ceiling piston complete the enclosure. Collisions are fully elastic
by default (`restitution = 1.0`).

### v1 — event-driven MD (exact)

- **Particle–parabola / walls**: contact time solved analytically as a
  quadratic; velocity reflected about the surface normal.
- **Particle–particle**: all pairwise contact times solved analytically
  (gravity cancels in the relative frame). The earliest event wins; all
  particles advance to that moment and exchange velocities via the elastic
  impulse formula.
- **Ceiling piston**: binary search for the contact time with the moving
  surface; the bounce reflects velocity relative to the piston.

### v2 — cell-list Verlet MD

Each outer timestep (1/120 s) is divided into `substeps` sub-steps (default 8):

1. **Verlet integrate** under gravity — O(N), OpenMP-parallel.
2. **Boundary conditions** — overshoot past a wall is *mirrored* back and the
   velocity reflected. Wall impulses (pressure) and piston work are
   accumulated here.
3. **Build cell list** — cell size = 2 × max radius, so overlapping pairs are
   always in the same or adjacent cells. Non-empty cells are bucketed into 9
   colours by (cx mod 3, cy mod 3).
4. **Resolve overlaps** — `passes` passes (default 4). Cells of one colour are
   ≥ 3 apart, so their 3×3 neighbourhoods are disjoint and every colour is
   processed in parallel with no locks. Each overlapping pair gets an elastic
   impulse (if approaching) and a positional push split by mass.

**Energy-consistent position corrections.** Any vertical move made outside
the integrator (wall mirror, overlap push, re-clamp) rescales the particle's
speed so that v² += 2·g·Δy — the move is energetically identical to free fall
over the same height. Without this, floor clamps under gravity inject energy
steadily (measured +1.74 % over 40 s at N = 1000, g = 500); with it, +0.002 %.
At g = 0 the drift is exactly zero, since impulses conserve KE exactly.

The initial lattice is stacked from the **floor up**, so the temperature is set
by `mono_speed` rather than by the drop (a 950 px fall at g = 500 would add
475 000 KE per unit mass), and the lattice sits clear of the piston at its
lowest point.

### Two-species gas

| Property | Light (cyan) | Heavy (orange) |
|----------|-------|-------|
| Mass | 1 | `mass_ratio` |
| Radius | `radius` | `radius × ∛mass_ratio` |

Collisions between unequal masses use `J = (1+e)·mᵢmⱼ/(mᵢ+mⱼ)·Δvₙ`.

### Equipartition, thermalisation, barometric profile

At equilibrium each species has the same mean KE per particle:
`½·m_light·⟨v²⟩ = ½·m_heavy·⟨v²⟩ = kT` (2D: mean KE per particle = kT).
Starting monoenergetic, the KE distribution relaxes to the 2D
Maxwell-Boltzmann exponential `f(KE) ∝ exp(−KE/kT)`. Under gravity each
species forms an atmosphere with scale height `kT/(m·g)`, so heavy particles
settle into a thinner layer.

### Oscillating piston

`y_ceil(t) = radius + ceil_amp × (1 + sin(2π × ceil_freq × t))`. For
near-zero heating the piston must be quasi-static:
`ceil_amp × 2π × ceil_freq ≪ mono_speed`. When faster, irreversible **Fermi
acceleration** heats the gas (measured: +137 % in 15 s at amp 30, 0.5 Hz,
speed 220 — ratio 0.43). The analysis panel shows the P–V loop and the work
done per cycle.

---

## Analysis panel (v2 only, right-hand column)

| Plot / readout | What it shows |
|---|---|
| **EQUIPARTITION KE/N** | Light and heavy mean KE per particle vs time; white line = overall mean. Press `M` and watch the two converge. |
| **H-FUNCTION** | `H = Σ p ln p` over fixed KE bins (KE marginal). Decreases toward the yellow equilibrium value `H_eq` computed from the exponential distribution at the current kT (Boltzmann H-theorem). |
| **DENSITY N(Y)** | Butterfly plot: light particles' vertical histogram grows left, heavy grows right. Yellow lines = barometric theory `n(y) ∝ exp(−m g h / kT) · width(y)`, normalised per species. Scale heights shown below. |
| **WALL PRESSURE** | `P = Σ|Δp| / (Δt · perimeter)` from wall impulses; `P_ideal = N kT / A`; `Z = PA/NkT` vs Henderson's hard-disc EOS `Z = (1 + η²/8)/(1−η)²`. Z is greyed out when gravity > 0 (the floor carries the weight). Measured at g = 0: Z = 1.03 vs 1.05 (N = 300), 1.60 vs 1.33 (N = 1000, η = 0.13). |
| **P–V LOOP** | Pressure vs vessel area, newest point white, trail fading. With the piston on, shows the cycle loop and the work per period (from piston–particle energy transfer, exact). |
| **HEAVY-PARTICLE MSD** | Ensemble mean-squared displacement of all heavy particles from their positions at the last reset (`T`). One heavy particle is ringed and trails on the main view. Ballistic → diffusive → saturated by the walls. |
| **KINETIC THEORY** | Collision rate per particle, mean speed, measured mean free path `⟨v⟩/ν`, and the Enskog prediction `λ = 1/(2√2 · n · σ · χ)` with `χ = (1 − 7η/16)/(1−η)²`, using the number-weighted local density from the profile bins. At g = 0: 18.6 vs 18.2 px (N = 1000). |
| **PERFORMANCE** | Physics ms per frame, FPS, thread count, substeps × passes. |

---

## Building

### Dependencies

- **SDL2** (≥ 2.0.18 for batched sprite rendering; older falls back automatically)
- **ffmpeg** — video export (optional)
- **OpenMP** — for v2 parallelism (optional; single-threaded without it)

```bash
# Linux
apt install libsdl2-dev ffmpeg

# macOS — Apple clang has no OpenMP; install Homebrew GCC (or libomp)
brew install sdl2 ffmpeg gcc
```

### Build

```bash
make          # builds ./ballsim (v1) and ./ballsim_v2 (v2)
make run2     # launch v2
make bench    # headless benchmark sweep N = 300 … 5000
make clean
```

The Makefile detects Homebrew `gcc-NN` on macOS and uses it for v2; it falls
back to clang + libomp, then to no OpenMP with a warning.

### WebAssembly (itch.io / browser) — v1 physics

```bash
source /path/to/emsdk/emsdk_env.sh
make -f Makefile.emscripten
python3 -m http.server 8080     # open http://localhost:8080/ballsim.html
mv ballsim.html index.html; mv ballsim.js index.js; mv ballsim.wasm index.wasm
zip ballsim.zip index.html index.js index.wasm
```

Upload to itch.io as an HTML project, embed size **980 × 960**, and enable
SharedArrayBuffer support.

---

## Command line (v2)

```
./ballsim_v2 [options]
  --config FILE          load config (default: ballsim_v2.cfg)
  --write-config [FILE]  write default config and exit
  --bench N STEPS        headless benchmark, no window
  --screenshot FILE K    run K steps offscreen, save BMP, exit
  --substeps K           Verlet sub-steps per 1/120 s (default 8)
  --passes P             overlap-resolution passes per sub-step (default 4)
  --threads T            OpenMP thread count
  --seed S               RNG seed
```

`substeps` and `passes` can also be set in the config file.

### Benchmark results

2-core container, `--bench N 600` (5 s of simulation), radius 4 unless noted:

| N | radius | substeps×passes | ms/step | realtime | drift / 5 s | max overlap |
|---|---|---|---|---|---|---|
| 300 | 4 | 8×4 | 0.9 | 9.1× | +0.0004 % | 0.5 % |
| 1000 | 4 | 8×4 | 2.1 | 4.0× | +0.002 % | 0.9 % |
| 2000 | 4 | 8×4 | 4.8 | 1.8× | +0.07 % | 0.8 % |
| 2000 | 4 | 4×4 | 2.5 | 3.4× | +0.07 % | 4.5 % |
| 5000 | 2 | 8×4 | 11.2 | 0.75× | +0.001 % | 3.5 % |
| 5000 | 2 | 4×4 | 6.0 | 1.4× | +0.004 % | 7.5 % |

"Realtime" is simulated seconds per wall second; ≥ 1 means 60 fps is
reachable. Drift is insensitive to the solver settings; only the residual
overlap (as % of touch distance — 4.5 % ≈ half a pixel) grows with fewer
substeps. On an 8-core machine expect roughly 2–3× these rates. Going
`--substeps 4` is the first thing to try when N is large.

---

## Configuration

```bash
./ballsim_v2 --write-config ballsim_v2.cfg
$EDITOR ballsim_v2.cfg
./ballsim_v2
```

| Parameter | v1 default | v2 default | Slider | Effect |
|-----------|-----------|-----------|--------|--------|
| `sim_w` | 700 | 700 | — | Simulation area width (restart required) |
| `sim_h` | 720 | 960 | — | Simulation area height (restart required) |
| `n_particles` | 150 | 300 | ✓ | Starting particle count (v1 ≤ 500, v2 ≤ 5000) |
| `gravity` | 500 | 500 | ✓ | px/s² |
| `radius` | 5 | 4 | ✓ | Base radius for light particles |
| `restitution` | 1.0 | 1.0 | ✓ | < 1.0 adds damping |
| `mono_speed` | 220 | 220 | ✓ | Speed for M-key reset; sets the temperature |
| `para_a` | 0.003 | 0.003 | ✓ | Parabola curvature; 0 = flat box |
| `ceil_amp` | 0 | 0 | ✓ | Piston amplitude; 0 = fixed ceiling |
| `ceil_freq` | 0.5 | 0.5 | ✓ | Piston frequency, Hz |
| `heavy_frac` | 0.5 | 0.5 | ✓ | Fraction of heavy particles |
| `mass_ratio` | 4.0 | 4.0 | ✓ | Heavy / light mass (1–20) |
| `substeps` | — | 8 | — | v2 solver sub-steps |
| `passes` | — | 4 | — | v2 resolution passes |

v2 window is `sim_w + 580` wide (two panels). Reduce `radius` as you raise
`n_particles`: N = 2000 at radius 3, N = 5000 at radius 2.

---

## Controls

| Key | Action |
|-----|--------|
| `M` | Reset — monoenergetic (equal KE per particle) |
| `R` | Reset — random speeds |
| `C` | Cycle colour mode: SPECIES → SPEED → KE → HEIGHT → RAINBOW (v2) |
| `T` | Reset the MSD reference and tracer trail (v2) |
| `V` | Toggle MP4 recording (Retina-aware in v2) |
| `+` / `=` | Add a particle (v2: spawned under the ceiling) |
| `-` | Remove a particle |
| `Q` / `Esc` | Quit |

Sliders: click a track to jump, drag to scrub. Changes affecting placement
trigger a reset.

---

## Stats panel (middle column)

| Readout | Meaning |
|---------|---------|
| TOTAL ENERGY (E/E0) | Energy normalised to the first sample — 1.0000 = conserved |
| KINETIC / POTENTIAL | Fractions of total |
| TEMPERATURE (KE/N) | Mean KE per particle = kT in 2D |
| LIGHT / HEAVY KE/N | Per-species mean KE — converge at equilibrium |
| COLLISIONS | Cumulative particle–particle collisions |
| ENERGY / E0 | Trace of the last 280 samples; green near 1.0, red when drifted |
| KE DISTRIBUTION vs MB | Measured histogram vs the 2D Maxwell-Boltzmann curve |

---

## Video recording

`V` starts/stops. Output `ballsim_v2_YYYYMMDD_HHMMSS.mp4`, H.264, 60 fps,
`-crf 18`, piped straight into ffmpeg; capture is gated on wall-clock time so
playback matches on-screen speed. Requires `ffmpeg` on `PATH`. Not available
in the WASM build.

---

## Repository layout

```
.
├── verlet_collisions_two_gasses_video.c   ← v1 source (event-driven, untouched)
├── ballsim_v2.c                           ← v2.1 source (cell-list MD + OpenMP + analysis)
├── ballsim_wasm.c                         ← WASM source (v1 physics, browser build)
├── Makefile                               ← builds v1 + v2; `make bench`
├── Makefile.emscripten                    ← builds WASM with emcc
├── shell.html                             ← HTML wrapper for the WASM build
└── README.md
```

---

## Suggested experiments

**Equipartition** — `M` with `mass_ratio = 8`. Watch the two lines in the
EQUIPARTITION plot converge on the white mean.

**H-theorem** — `M` (monoenergetic spike), watch H fall onto the yellow H_eq
line as the KE histogram relaxes to the exponential.

**Barometric atmosphere** — default settings. Heavy scale height is ¼ of the
light one; the yellow theory lines follow the measured profile. Raise
`gravity` to flatten both.

**Hard-disc equation of state** — `gravity = 0`, `para_a = 0`, single species
(`heavy_frac = 0`). Raise `n_particles` and compare Z with Henderson. Let the
lattice start settle for ~20 s first.

**Mean free path** — same setup; MFP MEASURED should track MFP ENSKOG within a
few per cent. Change radius and watch both move together.

**Brownian motion** — `mass_ratio = 20`, `heavy_frac = 0.05`, `gravity = 0`.
Press `T`; the MSD curve shows the ballistic-to-diffusive crossover before the
walls saturate it.

**Piston work** — `para_a = 0`, `gravity = 0`, `ceil_amp = 30`,
`ceil_freq = 0.5`. The P–V loop traces each cycle; W/CYC shows the energy
pumped in per period. Reduce amplitude × frequency until it approaches zero —
that is the quasi-static limit.

**Scale up** — drag PARTICLES toward 5000 with radius 2; `--substeps 4` if the
PHYSICS/FRAME readout exceeds ~14 ms.
