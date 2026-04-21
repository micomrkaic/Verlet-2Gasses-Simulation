# Parabolic Vessel — Hard Sphere Gas Simulation

A real-time 2D molecular dynamics simulation of hard-sphere particles in a
configurable vessel, written in C with SDL2. Demonstrates thermalisation,
energy conservation, Maxwell-Boltzmann statistics, equipartition between
species, and acoustic compression waves via an oscillating piston ceiling.
Supports video export to MP4.

---

## Physics

Particles are hard spheres under constant downward gravity inside a vessel
whose floor is a parabolic curve:

```
y_floor(x) = PARA_BASE - para_a * (x - sim_w/2)²
```

Setting `para_a = 0` gives a flat rectangular box. Vertical side walls and a
movable ceiling piston complete the enclosure.

**Collisions are fully elastic** by default (`restitution = 1.0`):

- **Particle–parabola / walls**: the contact time is solved analytically as a
  quadratic equation — no iteration needed. The trajectory under gravity is
  substituted into the parabola equation, yielding exact coefficients `A`, `B`,
  `C` whose roots give the precise contact moment. Velocity is then reflected
  about the surface normal.
- **Particle–particle**: fully event-driven MD — all pairwise contact times are
  solved analytically (gravity cancels in the relative frame, making the
  quadratic exact). The earliest event wins; all particles advance to that
  moment and exchange velocity components via the elastic impulse formula.
- **Ceiling piston**: binary search finds the exact contact time with the
  moving surface; the bounce reflects velocity relative to the piston,
  transferring the correct momentum.

Energy conservation is **< 0.15% drift over 12 000 frames** with a mixed
two-species gas, compared to ~80% drift from naive post-step overlap detection.

### Two-species gas

Particles come in two species — **light** (cyan) and **heavy** (orange). Their
properties are:

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

where `dvₙ` is the relative normal velocity at contact. This conserves both
momentum and kinetic energy exactly.

### Equipartition of energy

Starting from a monoenergetic state (all particles with equal KE), collisions
gradually drive the system toward thermal equilibrium. At equilibrium the
**equipartition theorem** requires each species to have the same mean KE per
particle:

```
½·m_light·⟨v_light²⟩  =  ½·m_heavy·⟨v_heavy²⟩  =  kT
```

Heavy particles slow down and light ones speed up until their mean KEs
converge. The LIGHT KE/N and HEAVY KE/N panel readouts let you watch this
in real time.

### Thermalisation

A single particle bouncing in the parabola conserves its own energy
indefinitely — the bowl is integrable. Thermalisation only occurs through
**particle–particle collisions**, which redistribute kinetic energy. Starting
monoenergetic (all same KE), the distribution broadens and converges to the
2D Maxwell-Boltzmann exponential:

```
f(KE) ∝ exp(−KE / kT)    where kT = mean KE per particle
```

### Oscillating piston

When `ceil_amp > 0` the ceiling acts as a sinusoidal piston:

```
y_ceil(t) = radius + ceil_amp × (1 + sin(2π × ceil_freq × t))
```

This drives compression waves that propagate visibly through the gas. For
near-zero net energy drift the piston must be **quasi-static**:

```
ceil_amp × 2π × ceil_freq  ≪  mono_speed   (ratio < ~0.1)
```

When the piston moves faster than the gas particles, irreversible **Fermi
acceleration** heats the gas permanently — this is real physics, not a bug.

---

## Building

### Dependencies

#### Required

| Dependency | Purpose | Minimum version |
|------------|---------|-----------------|
| **GCC or Clang** | C compiler | GCC 7 / Clang 6 (C11 support required) |
| **GNU Make** | Build system | Any recent version |
| **SDL2** | Window, renderer, input | 2.0.5+ |
| **libm** | Math library (`-lm`) | Bundled with glibc / system libc |

#### Optional

| Dependency | Purpose |
|------------|---------|
| **ffmpeg** | Video export — press `V` to record. Must be on your `PATH`. If absent, pressing `V` prints a warning and recording is silently skipped; the simulation continues normally. |
| **pkg-config** | Used by the Makefile to auto-detect SDL2 flags on Linux. If missing, the Makefile falls back to `-I/usr/include/SDL2 -lSDL2`. |

---

### Installing dependencies

#### Ubuntu / Debian

```bash
# Build tools
sudo apt install build-essential pkg-config

# SDL2 (runtime + development headers)
sudo apt install libsdl2-dev

# ffmpeg (optional — for video recording)
sudo apt install ffmpeg
```

#### Fedora / RHEL / CentOS

```bash
sudo dnf install gcc make pkgconf-pkg-config SDL2-devel

# ffmpeg (optional — may require RPM Fusion repo)
sudo dnf install ffmpeg
```

#### Arch Linux

```bash
sudo pacman -S base-devel sdl2

# ffmpeg (optional)
sudo pacman -S ffmpeg
```

#### macOS (Homebrew)

```bash
# Install Homebrew if needed: https://brew.sh
brew install sdl2

# ffmpeg (optional)
brew install ffmpeg
```

Xcode Command Line Tools provide `gcc` (Clang) and `make`:

```bash
xcode-select --install
```

---

### Building

A `Makefile` is included. From the project directory:

```bash
make          # compile → ./ballsim
make run      # compile and launch immediately
make clean    # remove the binary
```

The Makefile auto-detects the platform and SDL2 paths via `pkg-config` on
Linux or Homebrew paths on macOS. No manual flag editing is needed in the
normal case.

#### Manual build (if make is unavailable)

Linux:
```bash
gcc -O2 -std=c11 -Wall -Wextra -Wno-unused-parameter \
    $(pkg-config --cflags --libs sdl2) \
    -o ballsim verlet_collisions_two_gasses_video.c -lm
```

macOS (Apple Silicon):
```bash
gcc -O2 -std=c11 \
    -I/opt/homebrew/include/SDL2 -L/opt/homebrew/lib \
    -o ballsim verlet_collisions_two_gasses_video.c -lSDL2 -lm
```

Intel Mac: replace `/opt/homebrew` with `/usr/local`.

> **Include path note:** the source uses `#include <SDL2/SDL.h>`. If your
> SDL2 installation uses a flat include path, change it to `#include <SDL.h>`.

---

## Configuration

Parameters are loaded from a plain-text config file at startup. Most can
also be adjusted live via the **in-app sliders** in the right-hand panel.

### Quick start

```bash
./ballsim --write-config ballsim.cfg   # write a default template
$EDITOR ballsim.cfg                    # edit to taste
./ballsim                              # auto-loads ballsim.cfg
./ballsim --config experiment.cfg      # explicit path
```

### Box dimensions

`sim_w` and `sim_h` **must be set in the config file** — they require a window
resize and cannot be changed via sliders. Edit `ballsim.cfg` and restart.

```ini
sim_w = 900    # simulation area width  (200 – 2560 px)
sim_h = 800    # simulation area height (200 – 1440 px)
```

### Full parameter reference

| Parameter | Default | Slider | Effect |
|-----------|---------|--------|--------|
| `sim_w` | 900 | — | Simulation area width (restart required) |
| `sim_h` | 800 | — | Simulation area height (restart required) |
| `n_particles` | 150 | ✓ | Starting particle count (1–500) |
| `gravity` | 500 px/s² | ✓ | Stronger → faster bouncing, more PE at bottom |
| `radius` | 5 px | ✓ | Base radius for light particles; heavy scales up |
| `restitution` | 1.0 | ✓ | < 1.0 adds damping; energy trace slopes down |
| `mono_speed` | 220 px/s | ✓ | Reference speed for M-key reset (equal KE start) |
| `para_a` | 0.003 | ✓ | Parabola curvature; 0 = flat rectangular box |
| `ceil_amp` | 0 px | ✓ | Piston amplitude; 0 = fixed ceiling |
| `ceil_freq` | 0.5 Hz | ✓ | Piston oscillation frequency |
| `heavy_frac` | 0.5 | ✓ | Fraction of particles that are heavy (0–1) |
| `mass_ratio` | 4.0 | ✓ | Mass of heavy / mass of light (1–20) |

Parameters are clamped to safe ranges on load — typos won't crash the sim.
Unknown keys are silently ignored.

---

## In-app sliders

The lower portion of the right-hand panel contains ten sliders covering all
runtime-tunable parameters. Click anywhere on a track to jump to that value;
hold and drag to scrub continuously. The active slider highlights in bright
cyan. Changes that affect particle placement (radius, para_a, n_particles,
mono_speed, heavy_frac, mass_ratio) trigger an automatic reset.

---

## Controls

| Key | Action |
|-----|--------|
| `M` | Reset — monoenergetic (equal KE per particle, best for equipartition demo) |
| `R` | Reset — random speeds |
| `V` | Toggle video recording (starts/stops MP4 export) |
| `+` / `=` | Add a particle (up to 500) |
| `-` | Remove a particle |
| `Q` / `Esc` | Quit (also finalises any active recording) |

---

## Stats panel

### Numeric readouts

| Readout | Meaning |
|---------|---------|
| BOX SIZE | Current sim_w × sim_h in pixels |
| PARTICLES | Live particle count |
| TOTAL ENERGY (E/E0) | Energy normalised to initial — 1.0000 = perfectly conserved |
| KINETIC (KE/E) | Kinetic energy as fraction of total |
| POTENTIAL (PE/E) | Potential energy as fraction of total |
| TEMPERATURE (KE/N) | Mean KE per particle across all species |
| LIGHT KE/N | Mean KE per light particle — converges to HEAVY KE/N at equilibrium |
| HEAVY KE/N | Mean KE per heavy particle — converges to LIGHT KE/N at equilibrium |
| COLLISIONS | Cumulative particle–particle collision count |

KE/E and PE/E always sum to 1.0. At equilibrium in a parabolic bowl expect
KE/E ≈ 0.55–0.65 depending on how high the particles bounce.

Potential energy is computed as `PE = m·g·h` with the per-particle mass, so
both species contribute correctly to the total energy.

### REC indicator

Top-right corner of the panel shows recording status. When idle a dim grey dot
and `V-REC` hint are shown. When recording, a red dot blinks and the frame
counter increments. The indicator updates every display frame.

### Energy trace

Ring buffer of the last ~280 samples of E/E0. The white horizontal line marks
1.0 (perfect conservation). Bars colour from green (near 1.0) toward red as
drift grows. The DRIFT label shows the current percentage deviation.

### KE distribution

Histogram of per-particle kinetic energies (cyan bars, exponentially smoothed)
overlaid with the theoretical 2D Maxwell-Boltzmann curve (yellow). Press `M`
to reset to a monoenergetic spike and watch it relax to the exponential shape
over a few seconds as collisions redistribute energy.

### Species legend

Below the histogram: cyan swatch = light particles, orange swatch = heavy
particles.

---

## Video recording

Press `V` to start recording; press `V` again (or `Q`) to stop. The output
file is written to the current directory with a timestamped name:

```
ballsim_YYYYMMDD_HHMMSS.mp4
```

**Format:** MP4 with H.264 video (`-crf 18`, near-lossless), YUV420p pixel
format for broad compatibility (QuickTime, VLC, phones). Frame rate is 60 fps.

**Speed accuracy:** video capture is gated on wall-clock time (`SDL_GetTicks`),
not simulation time, so the recorded video plays back at exactly the same speed
you see on screen — independent of monitor refresh rate or CPU load.

**Requirements:** `ffmpeg` must be on your `PATH`. Frames are piped directly
into ffmpeg as raw RGB24 — no temporary files are written. If ffmpeg is not
installed or cannot be found, pressing `V` prints a warning to stderr and
recording is skipped — the simulation keeps running normally.

```bash
# Verify ffmpeg is available and on PATH
ffmpeg -version
```

---

## Suggested experiments

**Equipartition** — press `M` for a monoenergetic start with `mass_ratio = 8`
and `heavy_frac = 0.5`. Watch LIGHT KE/N and HEAVY KE/N in the panel: they
start equal, briefly diverge as lighter particles gain speed, then converge
back to the same value as collisions drive the system to equilibrium. The
convergence time shortens with higher particle count and larger radius.

**Mass segregation** — with high `mass_ratio` (try 16) and high `gravity`,
heavy particles preferentially settle near the bottom of the bowl while light
ones bounce higher — a gravitational analogue of atmospheric scale height.

**Thermalisation** — press `M` for monoenergetic start, watch the KE histogram
spike broaden into the exponential MB curve. The time scale depends on particle
density and radius.

**Flat box** — set `para_a = 0` for a rectangular box with uniform gravity.

**Compression waves** — set `para_a = 0`, `ceil_amp = 5`, `ceil_freq = 0.5`.
Watch density waves propagate downward from the piston on each downstroke.
Keep `ceil_amp × 2π × ceil_freq` well below `mono_speed` to avoid heating.

**Inelastic collapse** — set `restitution = 0.8` and watch kinetic energy
slowly drain into the energy trace. Particles cluster at the bottom.

**Dense gas** — increase `n_particles` and `radius` together; the collision
rate rises sharply and thermalisation becomes almost instantaneous.

**Record an experiment** — set up any of the above, press `M` to reset to a
clean monoenergetic state, press `V` to start recording, then watch the system
thermalise. Press `V` to stop. The resulting MP4 plays back at real-time speed.
