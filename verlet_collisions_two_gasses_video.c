#include <SDL2/SDL.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// ── Compile-time caps
// ─────────────────────────────────────────────────────────
#define MAX_PARTICLES 500
#define ENERGY_HIST 280 // samples in energy ring buffer
#define HIST_BINS 30    // bins in KE-distribution histogram
#define PANEL_W 280     // stats panel width — fixed

// ── Runtime config (all tunable parameters)
// ───────────────────────────────────
typedef struct {
  int n_particles;    // starting particle count
  double gravity;     // gravitational acceleration (px/s²)
  double radius;      // particle radius (px)
  double restitution; // coefficient of restitution (1.0 = elastic)
  double mono_speed;  // speed for monoenergetic reset (px/s)
  double para_a;      // parabola curvature coefficient
  int sim_w;          // simulation area width  (px)
  int sim_h;          // simulation area height (px)
  double ceil_amp;    // piston amplitude (px, 0 = fixed ceiling)
  double ceil_freq;   // piston oscillation frequency (Hz)
  double heavy_frac;  // fraction of particles that are heavy (0.0 – 1.0)
  double mass_ratio;  // mass of heavy / mass of light (>= 1.0)
} Config;

static Config cfg = {
    .n_particles = 150,
    .gravity = 500.0,
    .radius = 5.0,
    .restitution = 1.0,
    .mono_speed = 220.0,
    .para_a = 0.0030,
    .sim_w = 700,
    .sim_h = 720,
    .ceil_amp = 0.0,
    .ceil_freq = 0.5,
    .heavy_frac = 0.5,
    .mass_ratio = 4.0,
};

// Derived layout values (recomputed after config load)
static double CFG_PARA_BASE;
static int CFG_WIN_W;

static void config_derive(void) {
  // When curved, the parabola vertex sits 80px above the bottom so the bowl
  // rim doesn't clip the window edge.  When flat (para_a==0) the floor goes
  // all the way to the bottom of the sim area.
  CFG_PARA_BASE = (cfg.para_a > 1e-9) ? cfg.sim_h - 80.0 : cfg.sim_h;
  CFG_WIN_W = cfg.sim_w + PANEL_W;
}

// Parse a "key = value" config file.  Lines starting with '#' are comments.
// Unknown keys are silently ignored so old configs stay forward-compatible.
static void load_config(const char *path) {
  FILE *f = fopen(path, "r");
  if (!f) {
    fprintf(stderr, "Config: cannot open '%s', using defaults.\n", path);
    return;
  }
  char line[256];
  while (fgets(line, sizeof(line), f)) {
    // strip leading whitespace
    char *p = line;
    while (*p == ' ' || *p == '\t')
      p++;
    if (*p == '#' || *p == '\n' || *p == '\0')
      continue;
    char key[64];
    double val;
    if (sscanf(p, "%63s = %lf", key, &val) != 2)
      continue;
    if (!strcmp(key, "n_particles"))
      cfg.n_particles = (int)val;
    else if (!strcmp(key, "gravity"))
      cfg.gravity = val;
    else if (!strcmp(key, "radius"))
      cfg.radius = val;
    else if (!strcmp(key, "restitution"))
      cfg.restitution = val;
    else if (!strcmp(key, "mono_speed"))
      cfg.mono_speed = val;
    else if (!strcmp(key, "para_a"))
      cfg.para_a = val;
    else if (!strcmp(key, "sim_w"))
      cfg.sim_w = (int)val;
    else if (!strcmp(key, "sim_h"))
      cfg.sim_h = (int)val;
    else if (!strcmp(key, "ceil_amp"))
      cfg.ceil_amp = val;
    else if (!strcmp(key, "ceil_freq"))
      cfg.ceil_freq = val;
    else if (!strcmp(key, "heavy_frac"))
      cfg.heavy_frac = val;
    else if (!strcmp(key, "mass_ratio"))
      cfg.mass_ratio = val;
    else
      fprintf(stderr, "Config: unknown key '%s'\n", key);
  }
  fclose(f);
  // clamp to sane ranges
  if (cfg.n_particles < 1)
    cfg.n_particles = 1;
  if (cfg.n_particles > MAX_PARTICLES)
    cfg.n_particles = MAX_PARTICLES;
  if (cfg.gravity < 0)
    cfg.gravity = 0;
  if (cfg.radius < 1)
    cfg.radius = 1;
  if (cfg.radius > 20)
    cfg.radius = 20;
  if (cfg.restitution < 0)
    cfg.restitution = 0;
  if (cfg.restitution > 1)
    cfg.restitution = 1;
  if (cfg.mono_speed < 1)
    cfg.mono_speed = 1;
  if (cfg.para_a < 0)
    cfg.para_a = 0; // 0 = flat floor, allowed
  if (cfg.para_a > 0.01)
    cfg.para_a = 0.01;
  if (cfg.sim_w < 200)
    cfg.sim_w = 200;
  if (cfg.sim_w > 2560)
    cfg.sim_w = 2560;
  if (cfg.sim_h < 200)
    cfg.sim_h = 200;
  if (cfg.sim_h > 1440)
    cfg.sim_h = 1440;
  if (cfg.ceil_amp < 0)
    cfg.ceil_amp = 0;
  if (cfg.ceil_amp > 200)
    cfg.ceil_amp = 200;
  if (cfg.ceil_freq < 0)
    cfg.ceil_freq = 0;
  if (cfg.heavy_frac < 0) cfg.heavy_frac = 0;
  if (cfg.heavy_frac > 1) cfg.heavy_frac = 1;
  if (cfg.mass_ratio < 1) cfg.mass_ratio = 1;
  if (cfg.mass_ratio > 20) cfg.mass_ratio = 20;
}

// Write default config file so users have a template to start from
static void write_default_config(const char *path) {
  FILE *f = fopen(path, "w");
  if (!f) {
    fprintf(stderr, "Config: cannot write '%s'\n", path);
    return;
  }
  fprintf(f,
          "# Parabolic Vessel — Hard Sphere Gas  (config file)\n"
          "# All lines starting with '#' are comments.\n"
          "# Reload by restarting the simulation.\n"
          "\n"
          "sim_w       = %d       # simulation area width  px (200 – 2560)\n"
          "sim_h       = %d       # simulation area height px (200 – 1440)\n"
          "\n"
          "n_particles = %d       # starting count (1 – %d)\n"
          "gravity     = %.1f    # px/s²  (try 200 – 1000)\n"
          "radius      = %.1f     # px     (1 – 40)\n"
          "restitution = %.1f     # 0 = inelastic, 1 = elastic\n"
          "mono_speed  = %.1f   # px/s   initial speed for M-reset\n"
          "para_a      = %.4f   # parabola curvature (0 = flat floor, 0.001 – "
          "0.02 = curved)\n"
          "\n"
          "# Oscillating ceiling piston (set ceil_amp > 0 to enable)\n"
          "# Try: ceil_amp=60, ceil_freq=2.0\n"
          "ceil_amp    = %.1f     # piston amplitude px (0 = fixed ceiling)\n"
          "ceil_freq   = %.1f     # oscillation frequency Hz\n",
          cfg.sim_w, cfg.sim_h, cfg.n_particles, MAX_PARTICLES, cfg.gravity,
          cfg.radius, cfg.restitution, cfg.mono_speed, cfg.para_a, cfg.ceil_amp,
          cfg.ceil_freq);
  fclose(f);
  fprintf(stderr, "Config: wrote defaults to '%s'\n", path);
}

// ── Parabola vessel (defined before recording so capture_frame can use para_y)
static inline double para_y(double x) {
  double d = x - cfg.sim_w * 0.5;
  return CFG_PARA_BASE - cfg.para_a * d * d;
}
static inline double para_dydx(double x) {
  return -2.0 * cfg.para_a * (x - cfg.sim_w * 0.5);
}

// ── Oscillating ceiling (piston) ─────────────────────────────────────────────
static double sim_time = 0.0; // accumulated simulation time (seconds)

// y-coordinate of the ceiling at time t (uniform across full width = piston).
// The piston oscillates between cfg.radius (topmost) and cfg.radius +
// 2*ceil_amp. Centre of oscillation: cfg.radius + ceil_amp y_ceil(t) =
// cfg.radius + ceil_amp * (1 + sin(2π * ceil_freq * t))
//
// This means:
//   sin = -1  → y = cfg.radius          (piston at top, fully retracted)
//   sin =  0  → y = cfg.radius + amp    (piston at rest/midpoint)
//   sin = +1  → y = cfg.radius + 2*amp  (piston fully extended, maximum
//   compression)
//
// The piston NEVER rises above cfg.radius, so particles can never escape above
// it.
//
// QUASI-STATIC condition for near-zero net energy drift:
//   peak piston speed = ceil_amp * 2π * ceil_freq << typical particle speed
//   e.g. amp=5, freq=0.5Hz → peak=15.7 px/s vs 220 px/s (ratio 0.07) ✓
//   e.g. amp=60, freq=2.0Hz → peak=754 px/s vs 220 px/s (ratio 3.4) ✗ heats gas
static inline double ceil_y(double t) {
  if (cfg.ceil_amp < 1e-9)
    return cfg.radius;
  return cfg.radius +
         cfg.ceil_amp * (1.0 + sin(2.0 * M_PI * cfg.ceil_freq * t));
}

static inline double ceil_vy(double t) {
  if (cfg.ceil_amp < 1e-9)
    return 0.0;
  return cfg.ceil_amp * (2.0 * M_PI * cfg.ceil_freq) *
         cos(2.0 * M_PI * cfg.ceil_freq * t);
}

// ── Video recording ──────────────────────────────────────────────────────────
// Frames are piped as raw RGB24 into ffmpeg which encodes to H.264 MP4.
// ffmpeg must be installed: `apt install ffmpeg` / `brew install ffmpeg`.
// Output: ballsim_YYYYMMDD_HHMMSS.mp4  — one video frame per 1/REC_FPS
// seconds of simulation time, so playback speed always matches physics speed.
#define REC_FPS 60

static FILE  *rec_pipe     = NULL;
static bool   rec_active   = false;
static Uint32 rec_blink    = 0;
static long   rec_frames   = 0;
static Uint8 *rec_pixels   = NULL;
static char   rec_filename[64];
static Uint32 rec_wall_next = 0;    // wall-clock ms at which to emit next frame

static bool record_start(int win_w, int win_h) {
  if (rec_active) return true;
  time_t now = time(NULL);
  struct tm *lt = localtime(&now);
  strftime(rec_filename, sizeof(rec_filename),
           "ballsim_%Y%m%d_%H%M%S.mp4", lt);
  char cmd[512];
  snprintf(cmd, sizeof(cmd),
    "ffmpeg -y -f rawvideo -pixel_format rgb24"
    " -video_size %dx%d -framerate %d -i pipe:0"
    " -c:v libx264 -preset fast -crf 18"
    " -pix_fmt yuv420p"
    " \"%s\" 2>/dev/null",
    win_w, win_h, REC_FPS, rec_filename);
  rec_pipe = popen(cmd, "w");
  if (!rec_pipe) {
    fprintf(stderr, "Record: failed — is ffmpeg installed?\n");
    return false;
  }
  free(rec_pixels);
  rec_pixels = (Uint8 *)malloc((size_t)win_w * win_h * 3);
  if (!rec_pixels) { pclose(rec_pipe); rec_pipe = NULL; return false; }
  rec_active     = true;
  rec_frames     = 0;
  rec_wall_next  = SDL_GetTicks();
  fprintf(stderr, "Record: started → %s\n", rec_filename);
  return true;
}

static void record_stop(void) {
  if (!rec_active) return;
  rec_active = false;
  if (rec_pipe) { pclose(rec_pipe); rec_pipe = NULL; }
  fprintf(stderr, "Record: stopped — %ld frames → %s\n",
          rec_frames, rec_filename);
}

// Emit one video frame every 1/REC_FPS seconds of wall-clock time.
// Wall-clock gating (not sim-time) ensures the video plays back at
// exactly the same speed the user sees on screen.
static void capture_frame(SDL_Renderer *rend, int win_w, int win_h) {
  if (!rec_active || !rec_pipe || !rec_pixels) return;
  Uint32 now = SDL_GetTicks();
  if (now < rec_wall_next) return;
  SDL_RenderReadPixels(rend, NULL, SDL_PIXELFORMAT_RGB24,
                       rec_pixels, win_w * 3);
  fwrite(rec_pixels, 1, (size_t)win_w * win_h * 3, rec_pipe);
  rec_frames++;
  // Advance by exactly one frame interval, not by (now - rec_wall_next),
  // so jitter doesn't accumulate — any late frames shift the whole schedule.
  rec_wall_next += (Uint32)(1000 / REC_FPS);
  // If we've fallen more than one frame behind (heavy frame), re-sync
  // rather than trying to catch up with a burst of duplicate frames.
  if (now > rec_wall_next + (Uint32)(1000 / REC_FPS))
    rec_wall_next = now;
}

typedef struct {
  double x, y, vx, vy;
  double mass;    // 1.0 = light, mass_ratio = heavy
  double radius;  // per-particle radius (heavy gets larger)
  int species;    // 0 = light, 1 = heavy
  Uint8 r, g, b;
} Particle;
static Particle p[MAX_PARTICLES];

static double eng_buf[ENERGY_HIST]; // stores E/E0 (normalised)
static int eng_head = 0, eng_count = 0;
static double eng_max = 1.0, eng_min = 0.0;
static double eng_e0 = 0.0; // initial energy for normalisation

// Smoothed histogram (exponential moving average for visual stability)
static double hist_smooth[HIST_BINS];
static double hist_mb[HIST_BINS]; // theoretical Maxwell-Boltzmann curve
static double hist_ke_max = 1.0;  // current KE range for histogram x-axis

static long collision_count = 0; // total particle-particle collisions
// Per-species KE for stats panel
static double species_ke[2] = {0, 0};
static int    species_n[2]  = {0, 0};

static double randf(void) { return (double)rand() / ((double)RAND_MAX + 1.0); }

// ── Maxwell-Boltzmann distribution for 2D ideal gas ──────────────────────────
// In 2D, the speed distribution is f(v) = (m/kT)*v*exp(-mv²/2kT)
// In terms of KE = ½mv², f(KE) = (1/kT)*exp(-KE/kT)  (exponential)
// We normalise to [0,1] for display.
static void compute_mb_curve(double mean_ke) {
  if (mean_ke < 1.0)
    mean_ke = 1.0;
  double kT = mean_ke; // mean KE per particle = kT in 2D
  // Fill bins: each bin centre is at ke = (b+0.5)/HIST_BINS * hist_ke_max
  double peak = 0.0;
  double raw[HIST_BINS];
  for (int b = 0; b < HIST_BINS; b++) {
    double ke = (b + 0.5) / HIST_BINS * hist_ke_max;
    raw[b] = exp(-ke / kT);
    if (raw[b] > peak)
      peak = raw[b];
  }
  for (int b = 0; b < HIST_BINS; b++)
    hist_mb[b] = (peak > 0) ? raw[b] / peak : 0.0;
}

// ── Particle initialisation
// ───────────────────────────────────────────────────
static void init_particles(int n, bool monoenergetic) {
  // ── Pass 1: assign species, mass, radius, velocity, color ────────────────
  double max_r = cfg.radius; // track largest radius for placement
  for (int i = 0; i < n; i++) {
    int is_heavy = (randf() < cfg.heavy_frac) ? 1 : 0;
    p[i].species = is_heavy;
    p[i].mass    = is_heavy ? cfg.mass_ratio : 1.0;
    p[i].radius  = cfg.radius * cbrt(p[i].mass);
    if (p[i].radius > max_r) max_r = p[i].radius;

    double angle = randf() * 2.0 * M_PI;
    double speed = monoenergetic ? cfg.mono_speed : 80.0 + randf() * 300.0;
    // Equal KE start: ½mv² = ½·mono_speed²  →  v = mono_speed/√m
    double v = speed / sqrt(p[i].mass);
    p[i].vx = v * cos(angle);
    p[i].vy = v * sin(angle);

    double br = 0.75 + 0.25 * randf();
    if (is_heavy) {
      p[i].r = (Uint8)(255 * br);
      p[i].g = (Uint8)(110 * br);
      p[i].b = (Uint8)(20  * br);
    } else {
      p[i].r = (Uint8)(20  * br);
      p[i].g = (Uint8)(190 * br);
      p[i].b = (Uint8)(255 * br);
    }
  }

  // ── Pass 2: place on grid sized for the largest particle ─────────────────
  // Cell size = 2*max_r + gap so no two grid-adjacent particles ever overlap.
  double gap       = 2.0;                    // px breathing room between cells
  double cell      = 2.0 * max_r + gap;
  double piston_y0 = cfg.radius + cfg.ceil_amp; // piston midpoint at t=0
  double y_start   = piston_y0 + max_r + gap;   // first row top clearance
  double x_lo      = max_r + gap;
  double x_hi      = cfg.sim_w - max_r - gap;
  int cols = (int)((x_hi - x_lo) / cell);
  if (cols < 1) cols = 1;

  for (int i = 0; i < n; i++) {
    int col = i % cols, row = i / cols;
    p[i].x = x_lo + col * cell + max_r;      // centre of cell
    p[i].y = y_start + row * cell;
    // Tiny random jitter (well inside the cell) to break lattice symmetry
    p[i].x += (randf() - 0.5) * gap;
    p[i].y += (randf() - 0.5) * gap;
    // Hard clamp to sim area
    double ri = p[i].radius;
    if (p[i].x < ri + 1) p[i].x = ri + 1;
    if (p[i].x > cfg.sim_w - ri - 1) p[i].x = cfg.sim_w - ri - 1;
    if (p[i].y < ri + 1) p[i].y = ri + 1;
    if (p[i].y > cfg.sim_h - ri - 1) p[i].y = cfg.sim_h - ri - 1;
  }
  eng_head = eng_count = 0;
  eng_max = 1.0;
  eng_min = 0.0;
  eng_e0 = 0.0;
  for (int b = 0; b < HIST_BINS; b++) {
    hist_smooth[b] = 0.0;
    hist_mb[b] = 0.0;
  }
  hist_ke_max = 1.0;
  collision_count = 0;
  sim_time = 0.0;
}

// ── Energy measurement
// ────────────────────────────────────────────────────────
static double compute_energy(int n, double *ke_out, double *pe_out) {
  double ke = 0, pe = 0;
  species_ke[0] = species_ke[1] = 0;
  species_n[0]  = species_n[1]  = 0;
  for (int i = 0; i < n; i++) {
    double pke = 0.5 * p[i].mass * (p[i].vx * p[i].vx + p[i].vy * p[i].vy);
    ke += pke;
    pe += p[i].mass * cfg.gravity * (cfg.sim_h - p[i].radius - p[i].y);
    species_ke[p[i].species] += pke;
    species_n[p[i].species]++;
  }
  if (ke_out) *ke_out = ke;
  if (pe_out) *pe_out = pe;
  return ke + pe;
}

static void push_energy(double e) {
  if (eng_e0 == 0.0)
    eng_e0 = (e != 0.0) ? e : 1.0; // capture baseline on first push
  double norm = (eng_e0 != 0.0) ? e / eng_e0 : 1.0;
  eng_buf[eng_head] = norm;
  eng_head = (eng_head + 1) % ENERGY_HIST;
  if (eng_count < ENERGY_HIST)
    eng_count++;
  eng_max = eng_min = eng_buf[0];
  for (int i = 0; i < eng_count; i++) {
    if (eng_buf[i] > eng_max)
      eng_max = eng_buf[i];
    if (eng_buf[i] < eng_min)
      eng_min = eng_buf[i];
  }
  // Always include 1.0 in range; keep at least ±0.001 window around it
  if (eng_max < 1.0)
    eng_max = 1.0;
  if (eng_min > 1.0)
    eng_min = 1.0;
  if (eng_max - eng_min < 0.002) {
    eng_max = 1.0 + 0.001;
    eng_min = 1.0 - 0.001;
  }
}

// ── KE histogram update
// ───────────────────────────────────────────────────────
static void update_histogram(int n) {
  // Find max KE across all particles for axis scaling
  double max_ke = 1.0;
  for (int i = 0; i < n; i++) {
    double ke = 0.5 * p[i].mass * (p[i].vx * p[i].vx + p[i].vy * p[i].vy);
    if (ke > max_ke)
      max_ke = ke;
  }
  // Smooth the axis scale so it doesn't jump
  hist_ke_max += (max_ke * 1.1 - hist_ke_max) * 0.05;

  // Build raw histogram
  double raw[HIST_BINS] = {0};
  for (int i = 0; i < n; i++) {
    double ke = 0.5 * p[i].mass * (p[i].vx * p[i].vx + p[i].vy * p[i].vy);
    int b = (int)(ke / hist_ke_max * HIST_BINS);
    if (b >= HIST_BINS)
      b = HIST_BINS - 1;
    raw[b] += 1.0;
  }
  // Normalise raw to [0,1]
  double peak = 0;
  for (int b = 0; b < HIST_BINS; b++)
    if (raw[b] > peak)
      peak = raw[b];
  if (peak > 0)
    for (int b = 0; b < HIST_BINS; b++)
      raw[b] /= peak;
  // Exponential smoothing (α=0.15 → visual stability)
  for (int b = 0; b < HIST_BINS; b++)
    hist_smooth[b] = hist_smooth[b] * 0.85 + raw[b] * 0.15;

  // Recompute MB theoretical curve at current mean KE
  double ke_total = 0;
  for (int i = 0; i < n; i++)
    ke_total += 0.5 * p[i].mass * (p[i].vx * p[i].vx + p[i].vy * p[i].vy);
  compute_mb_curve(n > 0 ? ke_total / n : 1.0);
}

// ── CCD: advance one particle under gravity with parabola bounce
// ──────────────
static void advance_particle(Particle *q, double dt) {
  double ri = q->radius;
  double t_rem = dt;
  int iters = 0;
  while (t_rem > 1e-12 && iters++ < 8) {
    // Try free flight for full t_rem
    double tx  = q->x + q->vx * t_rem;
    double ty  = q->y + q->vy * t_rem + 0.5 * cfg.gravity * t_rem * t_rem;
    double tvy = q->vy + cfg.gravity * t_rem;
    if (ty + ri <= para_y(tx)) {
      q->x = tx; q->y = ty; q->vy = tvy;
      return;
    }
    // Analytic contact time: solve the quadratic
    //   (0.5g + para_a*vx²)·t² + (vy + 2·para_a·dx0·vx)·t + (y+r - PARA_BASE + para_a·dx0²) = 0
    // where dx0 = x - W/2.  Exact because gravity cancels cleanly with the
    // parabola's x²-dependence — no iteration needed.
    double dx0 = q->x - cfg.sim_w * 0.5;
    double A   = 0.5 * cfg.gravity + cfg.para_a * q->vx * q->vx;
    double B   = q->vy + 2.0 * cfg.para_a * dx0 * q->vx;
    double C   = q->y + ri - CFG_PARA_BASE + cfg.para_a * dx0 * dx0;
    double disc = B*B - 4.0*A*C;
    double tc;
    if (disc < 0) {
      tc = 0.0;   // inside surface, no real crossing — clamp
    } else {
      double sq = sqrt(disc);
      double t1 = (-B - sq) / (2.0*A);
      double t2 = (-B + sq) / (2.0*A);
      tc = t_rem + 1.0;
      if (t1 >= -1e-10 && t1 <= t_rem) tc = t1 < 0 ? 0 : t1;
      if (t2 >= -1e-10 && t2 <= t_rem && t2 < tc) tc = t2 < 0 ? 0 : t2;
      if (tc > t_rem) tc = 0.0;   // no root in range — clamp
    }
    // Advance to contact
    q->x  += q->vx * tc;
    q->y  += q->vy * tc + 0.5 * cfg.gravity * tc * tc;
    q->vy += cfg.gravity * tc;
    t_rem -= tc;
    // Reflect velocity about surface normal
    double dydx = para_dydx(q->x);
    double nx = -dydx, ny = 1.0, len = sqrt(nx*nx + ny*ny);
    nx /= len; ny /= len;
    double vn = q->vx * nx + q->vy * ny;
    if (vn > 0) {
      q->vx = (q->vx - 2.0 * vn * nx) * cfg.restitution;
      q->vy = (q->vy - 2.0 * vn * ny) * cfg.restitution;
    }
    // Push out by ri so next iteration's free-flight test starts outside
    q->y = para_y(q->x) - ri - ri * 1e-4;
  }
  // Consume any remaining time with free flight
  if (t_rem > 1e-12) {
    q->x  += q->vx  * t_rem;
    q->y  += q->vy  * t_rem + 0.5 * cfg.gravity * t_rem * t_rem;
    q->vy += cfg.gravity * t_rem;
  }
}

// ── Event-driven collision times
// ────────────────────────────────────────────── Screen coordinates: y
// increases DOWNWARD. Ceiling is near top (small y). Particle above ceiling
// means p.y < ceil_y. "gap" = p.y - ceil_y:  positive → particle below ceiling
// (normal)
//                         negative → particle above ceiling (inside piston)
// Particle approaches ceiling when its y is decreasing (vy < 0) faster than
// ceiling is decreasing (ceil_vy < 0 when piston rises).
// Relative velocity of particle w.r.t ceiling: vrel = p.vy - ceil_vy
//   vrel < 0 → particle moving toward ceiling (approaching)
//   vrel > 0 → particle moving away from ceiling (separating)
static double ceiling_time(int i, double dt) {
  double ri = p[i].radius;
  if (cfg.ceil_amp < 1e-9) {
    // Flat ceiling: particle top = p.y - ri must stay >= ri (ceiling at y=ri... actually ceiling at y=ri means particle center at y=2*ri)
    // ceiling surface is at y = ri (one base-radius from top), particle hits when p[i].y - ri = ri → p[i].y = 2*ri? No.
    // The ceiling is a hard wall at y = cfg.radius (one BASE radius from top, fixed).
    // Contact: particle center y = ceil_y + p[i].radius
    double ceil_surface = cfg.radius; // fixed ceiling y position
    double contact_y = ceil_surface + ri;
    if (p[i].vy >= 0) return dt + 1.0; // moving down, can't hit ceiling
    double a = 0.5 * cfg.gravity, b = p[i].vy, c = p[i].y - contact_y;
    double disc = b * b - 4.0 * a * c;
    if (disc < 0) return dt + 1.0;
    double sq = sqrt(disc);
    double t1 = (-b - sq) / (2.0 * a), t2 = (-b + sq) / (2.0 * a);
    double best = dt + 1.0;
    if (t1 >= -1e-10 && t1 <= dt) best = t1 < 0 ? 0 : t1;
    if (t2 >= -1e-10 && t2 <= dt && t2 < best) best = t2 < 0 ? 0 : t2;
    return best;
  }

  // Moving piston: contact when p[i].y = ceil_y(t) + ri
  double gap0 = p[i].y - (ceil_y(sim_time) + ri);
  if (gap0 < 0) return dt + 1.0; // already inside — clamp handles it
  double vrel = p[i].vy - ceil_vy(sim_time);
  if (vrel >= 0) return dt + 1.0;

  double ydt = p[i].y + p[i].vy * dt + 0.5 * cfg.gravity * dt * dt;
  double gapdt = ydt - (ceil_y(sim_time + dt) + ri);
  if (gapdt > 0) return dt + 1.0;

  double lo = 0, hi = dt;
  for (int k = 0; k < 52; k++) {
    double m = 0.5 * (lo + hi);
    double ym = p[i].y + p[i].vy * m + 0.5 * cfg.gravity * m * m;
    if (ym - (ceil_y(sim_time + m) + ri) > 0) lo = m;
    else hi = m;
  }
  return lo;
}

// Time until particle i hits the flat floor (y == CFG_PARA_BASE - cfg.radius)
// Only used when para_a == 0 (flat bottom); otherwise parabola CCD handles it.
// Solves: y + vy*t + 0.5*g*t^2 = floor_y
static double floor_time(int i, double dt) {
  if (cfg.para_a > 1e-9)
    return dt + 1.0;
  double floor_y = CFG_PARA_BASE - p[i].radius;
  if (p[i].vy <= 0) return dt + 1.0;
  double a = 0.5 * cfg.gravity, b = p[i].vy, c = p[i].y - floor_y;
  double disc = b * b - 4.0 * a * c;
  if (disc < 0)
    return dt + 1.0;
  double sq = sqrt(disc);
  double t1 = (-b - sq) / (2.0 * a), t2 = (-b + sq) / (2.0 * a);
  double best = dt + 1.0;
  if (t1 >= -1e-10 && t1 <= dt)
    best = t1 < 0 ? 0 : t1;
  if (t2 >= -1e-10 && t2 <= dt && t2 < best)
    best = t2 < 0 ? 0 : t2;
  return best;
}

static double lwall_time(int i, double dt) {
  if (p[i].vx >= 0) return dt + 1.0;
  double t = (p[i].radius - p[i].x) / p[i].vx;
  if (t >= -1e-10 && t <= dt) return t < 0 ? 0 : t;
  return dt + 1.0;
}

static double rwall_time(int i, double dt) {
  if (p[i].vx <= 0) return dt + 1.0;
  double t = (cfg.sim_w - p[i].radius - p[i].x) / p[i].vx;
  if (t >= -1e-10 && t <= dt) return t < 0 ? 0 : t;
  return dt + 1.0;
}

// Time until particles i,j first contact (gravity cancels in relative frame —
// exact)
static double contact_time(int i, int j, double dt) {
  double dx = p[i].x - p[j].x, dy = p[i].y - p[j].y;
  double dvx = p[i].vx - p[j].vx, dvy = p[i].vy - p[j].vy;
  double a = dvx * dvx + dvy * dvy;
  if (a < 1e-14)
    return dt + 1.0;
  double b = 2.0 * (dx * dvx + dy * dvy);
  double touch = p[i].radius + p[j].radius; // correct touch distance
  double c = dx * dx + dy * dy - touch * touch;
  if (c < 0) {
    return (b < 0) ? 0.0 : dt + 1.0;
  }
  double disc = b * b - 4.0 * a * c;
  if (disc < 0)
    return dt + 1.0;
  double t = (-b - sqrt(disc)) / (2.0 * a);
  if (t >= -1e-10 && t <= dt)
    return t < 0 ? 0 : t;
  return dt + 1.0;
}

// ── Full simulation step: event-driven ───────────────────────────────────────
// All events (ceiling, walls, particle-particle) compete for earliest time.
// Advance all particles to that time, process the event, repeat.
// This guarantees exact collision positions regardless of timestep size.
#define EV_PP 0
#define EV_CEIL 1
#define EV_LWALL 2
#define EV_RWALL 3
#define EV_FLOOR 4

static void step_sim(int n, double dt) {
  double t_done = 0;
  int max_events = n * 12;

  while (t_done < dt - 1e-12 && max_events-- > 0) {
    double t_rem = dt - t_done;
    double t_min = t_rem;
    int ev_type = -1, ei = -1, ej = -1;

    // Ceiling, floor (flat only), and wall events
    for (int i = 0; i < n; i++) {
      double tc = ceiling_time(i, t_rem);
      if (tc < t_min) {
        t_min = tc;
        ev_type = EV_CEIL;
        ei = i;
      }
      tc = floor_time(i, t_rem);
      if (tc < t_min) {
        t_min = tc;
        ev_type = EV_FLOOR;
        ei = i;
      }
      if (p[i].y < para_y(p[i].radius)) { // side walls only above parabola rim
        tc = lwall_time(i, t_rem);
        if (tc < t_min) {
          t_min = tc;
          ev_type = EV_LWALL;
          ei = i;
        }
        tc = rwall_time(i, t_rem);
        if (tc < t_min) {
          t_min = tc;
          ev_type = EV_RWALL;
          ei = i;
        }
      }
    }
    // Particle-particle contact events
    for (int i = 0; i < n - 1; i++) {
      for (int j = i + 1; j < n; j++) {
        double tc = contact_time(i, j, t_rem);
        if (tc < t_min) {
          t_min = tc;
          ev_type = EV_PP;
          ei = i;
          ej = j;
        }
      }
    }

    // Advance all particles to earliest event time
    for (int i = 0; i < n; i++) {
      advance_particle(&p[i], t_min);
      // Hard clamp: if particle ended up above ceiling (can happen when
      // advance_particle doesn't know about the moving piston), push it
      // back to the ceiling surface and reflect velocity.
      double cy = ceil_y(sim_time + t_done);
      if (p[i].y < cy + p[i].radius) {
        double vcy = ceil_vy(sim_time + t_done);
        double rel = p[i].vy - vcy;
        if (rel < 0)
          p[i].vy = -rel * cfg.restitution + vcy;
        p[i].y = cy + p[i].radius + 1e-4;
      }
    }
    t_done += t_min;

    // Process event — sim_time+t_done is the correct absolute time of this
    // event
    if (ev_type == EV_CEIL) {
      double ct = sim_time + t_done;
      double vcy = ceil_vy(ct);
      double rel_vy = p[ei].vy - vcy;
      if (rel_vy < 0) {
        p[ei].vy = (-rel_vy) * cfg.restitution + vcy;
      }
      p[ei].y = ceil_y(ct) + p[ei].radius + 1e-4;
    } else if (ev_type == EV_FLOOR) {
      p[ei].y = CFG_PARA_BASE - p[ei].radius;
      p[ei].vy = -fabs(p[ei].vy) * cfg.restitution;
    } else if (ev_type == EV_LWALL) {
      p[ei].x = p[ei].radius;
      p[ei].vx = fabs(p[ei].vx) * cfg.restitution;
    } else if (ev_type == EV_RWALL) {
      p[ei].x = cfg.sim_w - p[ei].radius;
      p[ei].vx = -fabs(p[ei].vx) * cfg.restitution;
    } else if (ev_type == EV_PP) {
      double dx = p[ej].x - p[ei].x, dy = p[ej].y - p[ei].y;
      double dist = sqrt(dx * dx + dy * dy);
      double touch = p[ei].radius + p[ej].radius;
      if (dist < 1e-12) dist = touch;
      double nx = dx / dist, ny = dy / dist;
      double dvn = (p[ej].vx - p[ei].vx) * nx + (p[ej].vy - p[ei].vy) * ny;
      if (dvn < 0) {
        double mi = p[ei].mass, mj = p[ej].mass;
        double J = 2.0 * mi * mj / (mi + mj) * dvn;
        p[ei].vx += (J / mi) * nx;
        p[ei].vy += (J / mi) * ny;
        p[ej].vx -= (J / mj) * nx;
        p[ej].vy -= (J / mj) * ny;
        collision_count++;
      }
      // push apart by one small step along normal to prevent re-collision
      double overlap = touch - dist;
      double push = (overlap > 0 ? overlap * 0.5 : 0) + 1e-4;
      p[ei].x -= push * nx;
      p[ei].y -= push * ny;
      p[ej].x += push * nx;
      p[ej].y += push * ny;
    }
  }
  sim_time += dt;
}

// ── Drawing
// ───────────────────────────────────────────────────────────────────
static void draw_particle(SDL_Renderer *rend, double x, double y, double rad,
                          Uint8 r, Uint8 g, Uint8 b) {
  int cx = (int)x, cy = (int)y, R = (int)rad;
  SDL_SetRenderDrawBlendMode(rend, SDL_BLENDMODE_BLEND);
  for (int dy2 = -R - 3; dy2 <= R + 3; dy2++)
    for (int dx2 = -R - 3; dx2 <= R + 3; dx2++) {
      double d = sqrt((double)(dx2 * dx2 + dy2 * dy2));
      if (d < R + 3 && d >= R) {
        double a = 1.0 - (d - R) / 3.0;
        SDL_SetRenderDrawColor(rend, r, g, b, (Uint8)(a * 70));
        SDL_RenderDrawPoint(rend, cx + dx2, cy + dy2);
      }
    }
  SDL_SetRenderDrawBlendMode(rend, SDL_BLENDMODE_NONE);
  for (int dy2 = -R; dy2 <= R; dy2++)
    for (int dx2 = -R; dx2 <= R; dx2++)
      if (dx2 * dx2 + dy2 * dy2 <= R * R) {
        double br = 1.0 - 0.3 * (dx2 * dx2 + dy2 * dy2) / (double)(R * R);
        SDL_SetRenderDrawColor(rend, (Uint8)(r * br), (Uint8)(g * br),
                               (Uint8)(b * br), 255);
        SDL_RenderDrawPoint(rend, cx + dx2, cy + dy2);
      }
}

static void draw_vessel(SDL_Renderer *rend) {
  SDL_SetRenderDrawColor(rend, 0x55, 0x77, 0xAA, 0xFF);
  // Floor: parabola curve (or flat line when para_a==0)
  for (int xs = 0; xs < cfg.sim_w - 1; xs++) {
    int y0 = (int)para_y(xs), y1 = (int)para_y(xs + 1);
    for (int t = -1; t <= 1; t++)
      SDL_RenderDrawLine(rend, xs, y0 + t, xs + 1, y1 + t);
  }
  // Ceiling piston: filled rectangle from top of screen down to piston face
  {
    int cy = (int)ceil_y(sim_time);
    SDL_SetRenderDrawColor(rend, 0x1A, 0x22, 0x3A,
                           0xFF); // dark fill (wall body)
    SDL_Rect piston = {0, 0, cfg.sim_w, cy};
    SDL_RenderFillRect(rend, &piston);
    // Bright edge at piston face
    SDL_SetRenderDrawColor(rend, 0x55, 0x77, 0xAA, 0xFF);
    for (int t = -1; t <= 1; t++)
      SDL_RenderDrawLine(rend, 0, cy + t, cfg.sim_w, cy + t);
  }
  // Side walls
  double rim_y = (cfg.para_a > 1e-9) ? para_y(cfg.radius) : CFG_PARA_BASE;
  SDL_Rect lw = {0, 0, 2, (int)rim_y}, rw = {cfg.sim_w - 2, 0, 2, (int)rim_y};
  SDL_RenderFillRect(rend, &lw);
  SDL_RenderFillRect(rend, &rw);
}

// ── Tiny 3×5 pixel font — digits 0-9, uppercase A-Z, '.', '-', '/', '(', ')'
// ── Each glyph is 5 rows of 3 bits (MSB = left column). Index: 0-9 = digits,
// 10-35 = A-Z, 36='.', 37='-', 38='/', 39='(', 40=')', 41=' '
static const Uint8 F[][5] = {
    // 0-9
    {0x7, 0x5, 0x5, 0x5, 0x7},
    {0x2, 0x6, 0x2, 0x2, 0x7},
    {0x7, 0x1, 0x7, 0x4, 0x7},
    {0x7, 0x1, 0x7, 0x1, 0x7},
    {0x5, 0x5, 0x7, 0x1, 0x1},
    {0x7, 0x4, 0x7, 0x1, 0x7},
    {0x7, 0x4, 0x7, 0x5, 0x7},
    {0x7, 0x1, 0x1, 0x1, 0x1},
    {0x7, 0x5, 0x7, 0x5, 0x7},
    {0x7, 0x5, 0x7, 0x1, 0x7},
    // A-Z (index 10-35)
    {0x2, 0x5, 0x7, 0x5, 0x5}, // A
    {0x6, 0x5, 0x6, 0x5, 0x6}, // B
    {0x7, 0x4, 0x4, 0x4, 0x7}, // C
    {0x6, 0x5, 0x5, 0x5, 0x6}, // D
    {0x7, 0x4, 0x6, 0x4, 0x7}, // E
    {0x7, 0x4, 0x6, 0x4, 0x4}, // F
    {0x7, 0x4, 0x5, 0x5, 0x7}, // G
    {0x5, 0x5, 0x7, 0x5, 0x5}, // H
    {0x7, 0x2, 0x2, 0x2, 0x7}, // I
    {0x1, 0x1, 0x1, 0x5, 0x7}, // J
    {0x5, 0x5, 0x6, 0x5, 0x5}, // K
    {0x4, 0x4, 0x4, 0x4, 0x7}, // L
    {0x5, 0x7, 0x7, 0x5, 0x5}, // M
    {0x5, 0x7, 0x5, 0x5, 0x5}, // N
    {0x2, 0x5, 0x5, 0x5, 0x2}, // O
    {0x6, 0x5, 0x6, 0x4, 0x4}, // P
    {0x2, 0x5, 0x5, 0x7, 0x3}, // Q
    {0x6, 0x5, 0x6, 0x5, 0x5}, // R
    {0x7, 0x4, 0x7, 0x1, 0x7}, // S
    {0x7, 0x2, 0x2, 0x2, 0x2}, // T
    {0x5, 0x5, 0x5, 0x5, 0x7}, // U
    {0x5, 0x5, 0x5, 0x5, 0x2}, // V
    {0x5, 0x5, 0x7, 0x7, 0x5}, // W
    {0x5, 0x5, 0x2, 0x5, 0x5}, // X
    {0x5, 0x5, 0x2, 0x2, 0x2}, // Y
    {0x7, 0x1, 0x2, 0x4, 0x7}, // Z
    // punctuation (index 36+)
    {0x0, 0x0, 0x0, 0x0, 0x2}, // 36 '.'
    {0x0, 0x0, 0x7, 0x0, 0x0}, // 37 '-'
    {0x1, 0x1, 0x2, 0x4, 0x4}, // 38 '/'
    {0x3, 0x4, 0x4, 0x4, 0x3}, // 39 '('
    {0x6, 0x1, 0x1, 0x1, 0x6}, // 40 ')'
    {0x0, 0x0, 0x0, 0x0, 0x0}, // 41 ' '
};
static int gl(char c) {
  if (c >= '0' && c <= '9')
    return c - '0';
  if (c >= 'A' && c <= 'Z')
    return 10 + (c - 'A');
  if (c >= 'a' && c <= 'z')
    return 10 + (c - 'a'); // lowercase maps to uppercase
  if (c == '.')
    return 36;
  if (c == '-')
    return 37;
  if (c == '/')
    return 38;
  if (c == '(')
    return 39;
  if (c == ')')
    return 40;
  return 41; // space / unknown
}
static void ds(SDL_Renderer *rend, const char *s, int x, int y, int sc, Uint8 r,
               Uint8 g, Uint8 b) {
  SDL_SetRenderDrawColor(rend, r, g, b, 255);
  for (int ci = 0; s[ci]; ci++)
    for (int row = 0; row < 5; row++)
      for (int col = 0; col < 3; col++)
        if (F[gl(s[ci])][row] & (0x4 >> col)) {
          SDL_Rect px = {x + ci * (3 * sc + sc) + col * sc, y + row * sc, sc,
                         sc};
          SDL_RenderFillRect(rend, &px);
        }
}

static int g_slider_panel_top = 500; // set by draw_panel each frame

// ── Slider system
// ───────────────────────────────────────────────────────────── Sliders can
// point to either a double or an int field in cfg.
typedef struct {
  const char *label;
  double min, max;
  double *dvalue; // non-NULL for double cfg fields
  int *ivalue;    // non-NULL for int cfg fields
  int needs_reset;
} Slider;

static Slider sliders[10];
static int n_sliders = 0;
static int active_slider = -1;

static void sliders_init(void) {
  n_sliders = 0;
#define SD(lbl, mn, mx, fld, rst)                                              \
  sliders[n_sliders++] = (Slider){lbl, mn, mx, &cfg.fld, NULL, rst}
#define SI(lbl, mn, mx, fld, rst)                                              \
  sliders[n_sliders++] = (Slider){lbl, mn, mx, NULL, &cfg.fld, rst}
  SI("PARTICLES", 1, 500, n_particles, 1);
  SD("GRAVITY", 0, 2000, gravity, 0);
  SD("RADIUS", 1, 20, radius, 1);
  SD("RESTITUTION", 0, 1, restitution, 0);
  SD("MONO SPEED", 10, 1000, mono_speed, 1);
  SD("PARA A", 0, 0.01, para_a, 1);
  SD("PISTON AMP", 0, 150, ceil_amp, 0);
  SD("PISTON FREQ", 0, 10, ceil_freq, 0);
  SD("HEAVY FRAC", 0, 1, heavy_frac, 1);
  SD("MASS RATIO", 1, 20, mass_ratio, 1);
#undef SD
#undef SI
}

static double slider_get(int i) {
  return sliders[i].dvalue ? *sliders[i].dvalue : (double)*sliders[i].ivalue;
}
static void slider_put(int i, double v) {
  if (sliders[i].dvalue)
    *sliders[i].dvalue = v;
  else
    *sliders[i].ivalue = (int)round(v);
}

static SDL_Rect slider_track(int i, int top) {
  // 18px header gap + 30px per row
  return (SDL_Rect){cfg.sim_w + 8, top + 18 + i * 30 + 16, PANEL_W - 16, 5};
}

static void draw_sliders(SDL_Renderer *rend, int top) {
  ds(rend, "PARAMETERS", cfg.sim_w + 8, top, 2, 0xAA, 0xBB, 0xCC);
  for (int i = 0; i < n_sliders; i++) {
    Slider *s = &sliders[i];
    SDL_Rect tr = slider_track(i, top);
    int ly = tr.y - 13;
    // Label at scale 2
    ds(rend, s->label, tr.x, ly, 2, 0x88, 0x99, 0xBB);
    // Value at scale 2, right-aligned (each char = 3*2+2 = 8px wide)
    char buf[16];
    double v = slider_get(i);
    if (s->ivalue)
      snprintf(buf, sizeof(buf), "%d", (int)v);
    else if (s->max <= 0.1)
      snprintf(buf, sizeof(buf), "%.4f", v);
    else if (s->max <= 1.0)
      snprintf(buf, sizeof(buf), "%.3f", v);
    else if (s->max <= 10.0)
      snprintf(buf, sizeof(buf), "%.2f", v);
    else
      snprintf(buf, sizeof(buf), "%.1f", v);
    ds(rend, buf, tr.x + tr.w - (int)strlen(buf) * 8, ly, 2, 0xFF, 0xDD, 0x88);
    // Track
    SDL_SetRenderDrawColor(rend, 0x22, 0x2A, 0x3A, 0xFF);
    SDL_RenderFillRect(rend, &tr);
    double frac = (v - s->min) / (s->max - s->min);
    if (frac < 0)
      frac = 0;
    if (frac > 1)
      frac = 1;
    SDL_Rect fill = tr;
    fill.w = (int)(frac * tr.w);
    SDL_SetRenderDrawColor(rend, i == active_slider ? 0x88 : 0x44,
                           i == active_slider ? 0xCC : 0x88,
                           i == active_slider ? 0xFF : 0xCC, 0xFF);
    if (fill.w > 0)
      SDL_RenderFillRect(rend, &fill);
    SDL_Rect thumb = {tr.x + fill.w - 4, tr.y - 4, 8, 13};
    SDL_SetRenderDrawColor(rend, 0xDD, 0xEE, 0xFF, 0xFF);
    SDL_RenderFillRect(rend, &thumb);
  }
}

static int slider_hit(int mx, int my, int top) {
  for (int i = 0; i < n_sliders; i++) {
    SDL_Rect tr = slider_track(i, top);
    SDL_Rect h = {tr.x, tr.y - 6, tr.w, tr.h + 12};
    if (mx >= h.x && mx < h.x + h.w && my >= h.y && my < h.y + h.h)
      return i;
  }
  return -1;
}

static int slider_set(int idx, int mx, int top) {
  Slider *s = &sliders[idx];
  SDL_Rect tr = slider_track(idx, top);
  double frac = (double)(mx - tr.x) / tr.w;
  if (frac < 0)
    frac = 0;
  if (frac > 1)
    frac = 1;
  slider_put(idx, s->min + frac * (s->max - s->min));
  config_derive();
  return s->needs_reset;
}

static int draw_stat(SDL_Renderer *rend, int px, int y, const char *label,
                     const char *val, Uint8 vr, Uint8 vg, Uint8 vb) {
  ds(rend, label, px, y, 2, 0xAA, 0xBB, 0xCC);
  ds(rend, val, px, y + 12, 3, vr, vg, vb);
  return y + 30;
}

// ── Stats panel
// ───────────────────────────────────────────────────────────────
static void draw_panel(SDL_Renderer *rend, int n, double total_e, double ke,
                       double pe) {
  SDL_SetRenderDrawColor(rend, 0x0C, 0x0F, 0x1C, 0xFF);
  SDL_Rect bg = {cfg.sim_w, 0, PANEL_W, cfg.sim_h};
  SDL_RenderFillRect(rend, &bg);
  SDL_SetRenderDrawColor(rend, 0x44, 0x55, 0x77, 0xFF);
  SDL_RenderDrawLine(rend, cfg.sim_w, 0, cfg.sim_w, cfg.sim_h);

  int px = cfg.sim_w + 8, ty = 10;
  char buf[32];

  // ── REC indicator (top-right of panel) ───────────────────────────────────
  {
    int rx = cfg.sim_w + PANEL_W - 8;   // right-align
    if (rec_active) {
      // Blink: on 600 ms, off 300 ms
      bool dot_on = ((SDL_GetTicks() - rec_blink) % 900) < 600;
      if (dot_on) {
        SDL_SetRenderDrawColor(rend, 0xEE, 0x22, 0x22, 0xFF);
        // Filled circle approximation (7px radius)
        for (int dy2 = -5; dy2 <= 5; dy2++)
          for (int dx2 = -5; dx2 <= 5; dx2++)
            if (dx2*dx2 + dy2*dy2 <= 25)
              SDL_RenderDrawPoint(rend, rx - 6 + dx2, ty + 5 + dy2);
      }
      ds(rend, "REC", rx - 38, ty, 2, 0xEE, 0x33, 0x33);
      // Frame counter below
      snprintf(buf, sizeof(buf), "%ldf", rec_frames);
      ds(rend, buf, rx - (int)strlen(buf)*8, ty + 14, 2, 0xAA, 0x44, 0x44);
    } else {
      // Idle: dim grey dot + "V-REC" hint
      SDL_SetRenderDrawColor(rend, 0x33, 0x33, 0x33, 0xFF);
      for (int dy2 = -4; dy2 <= 4; dy2++)
        for (int dx2 = -4; dx2 <= 4; dx2++)
          if (dx2*dx2 + dy2*dy2 <= 16)
            SDL_RenderDrawPoint(rend, rx - 6 + dx2, ty + 4 + dy2);
      ds(rend, "V-REC", rx - 46, ty, 2, 0x44, 0x44, 0x55);
    }
  }

  // Box dimensions — set in ballsim.cfg (sim_w / sim_h)
  snprintf(buf, sizeof(buf), "%dx%d", cfg.sim_w, cfg.sim_h);
  ty = draw_stat(rend, px, ty, "BOX SIZE (cfg sim-w/h)", buf, 0x88, 0x99, 0xAA);
  ty += 4;

  snprintf(buf, sizeof(buf), "%d", n);
  ty = draw_stat(rend, px, ty, "PARTICLES", buf, 0xFF, 0xFF, 0xFF);
  ty += 4;

  snprintf(buf, sizeof(buf), "%.4f", (eng_e0 > 0) ? total_e / eng_e0 : 1.0);
  ty = draw_stat(rend, px, ty, "TOTAL ENERGY (E/E0)", buf, 0xFF, 0xCC, 0x44);
  ty += 4;

  snprintf(buf, sizeof(buf), "%.4f", (total_e != 0) ? ke / total_e : 0.0);
  ty = draw_stat(rend, px, ty, "KINETIC (KE/E)", buf, 0x44, 0xFF, 0xAA);
  ty += 4;

  snprintf(buf, sizeof(buf), "%.4f", (total_e != 0) ? pe / total_e : 0.0);
  ty = draw_stat(rend, px, ty, "POTENTIAL (PE/E)", buf, 0x88, 0xAA, 0xFF);
  ty += 4;

  double mean_ke = (n > 0) ? ke / n : 0.0;
  snprintf(buf, sizeof(buf), "%.1f", mean_ke);
  ty = draw_stat(rend, px, ty, "TEMPERATURE (KE/N)", buf, 0xFF, 0x66, 0x44);
  ty += 4;

  // Per-species mean KE — shows equipartition approach
  {
    double lke = species_n[0] > 0 ? species_ke[0] / species_n[0] : 0.0;
    double hke = species_n[1] > 0 ? species_ke[1] / species_n[1] : 0.0;
    snprintf(buf, sizeof(buf), "%.1f", lke);
    ty = draw_stat(rend, px, ty, "LIGHT KE/N", buf, 0x33, 0xCC, 0xFF);
    ty += 4;
    snprintf(buf, sizeof(buf), "%.1f", hke);
    ty = draw_stat(rend, px, ty, "HEAVY KE/N", buf, 0xFF, 0x88, 0x33);
    ty += 4;
  }

  snprintf(buf, sizeof(buf), "%ld", collision_count);
  ty = draw_stat(rend, px, ty, "COLLISIONS", buf, 0xCC, 0x88, 0xFF);
  ty += 10;

  // ── Energy trace ──────────────────────────────────────────────────────────
  ds(rend, "ENERGY / E0", px, ty, 2, 0xAA, 0xBB, 0xCC);
  ty += 13;
  int gx = cfg.sim_w + 4, gw = PANEL_W - 8, gh = 80, gy = ty;
  SDL_SetRenderDrawColor(rend, 0x06, 0x08, 0x10, 0xFF);
  SDL_Rect gb = {gx, gy, gw, gh};
  SDL_RenderFillRect(rend, &gb);

  double range = eng_max - eng_min;
// y pixel for a normalised value v
#define TRACE_Y(v) (gy + gh - 1 - (int)(((v) - eng_min) / range * (gh - 2)))

  // Dim grid lines at 25% intervals
  SDL_SetRenderDrawBlendMode(rend, SDL_BLENDMODE_BLEND);
  SDL_SetRenderDrawColor(rend, 0x22, 0x33, 0x55, 80);
  for (int f = 1; f <= 3; f++)
    SDL_RenderDrawLine(rend, gx, gy + gh - gh * f / 4, gx + gw,
                       gy + gh - gh * f / 4);
  SDL_SetRenderDrawBlendMode(rend, SDL_BLENDMODE_NONE);

  // Reference line at E/E0 = 1.0
  int ref_y = TRACE_Y(1.0);
  if (ref_y >= gy && ref_y < gy + gh) {
    SDL_SetRenderDrawBlendMode(rend, SDL_BLENDMODE_BLEND);
    SDL_SetRenderDrawColor(rend, 0xFF, 0xFF, 0xFF, 60);
    SDL_RenderDrawLine(rend, gx, ref_y, gx + gw, ref_y);
    SDL_SetRenderDrawBlendMode(rend, SDL_BLENDMODE_NONE);
  }

  // Trace bars
  if (eng_count > 1) {
    int cols = (eng_count < gw) ? eng_count : gw;
    for (int col = 0; col < cols; col++) {
      int idx = (eng_head - 1 - col + ENERGY_HIST) % ENERGY_HIST;
      double v = eng_buf[idx];
      double frac = (range > 0) ? (v - eng_min) / range : 0.5;
      int bar_y = TRACE_Y(v);
      if (bar_y < gy)
        bar_y = gy;
      if (bar_y >= gy + gh)
        bar_y = gy + gh - 1;
      int bx = gx + gw - 1 - col;
      // colour: green when near 1.0, red when drifted
      double dev = fabs(v - 1.0) / (range * 0.5 + 1e-9);
      if (dev > 1.0)
        dev = 1.0;
      Uint8 cr = (Uint8)(dev * 255);
      Uint8 cg = (Uint8)((1.0 - dev) * 200);
      SDL_SetRenderDrawColor(rend, cr, cg, 0x18, 0xFF);
      SDL_RenderDrawLine(rend, bx, bar_y, bx, gy + gh - 1);
      (void)frac;
    }
  }
  SDL_SetRenderDrawColor(rend, 0x2A, 0x3A, 0x55, 0xFF);
  SDL_RenderDrawRect(rend, &gb);

  // Current deviation label below graph
  if (eng_count > 0) {
    double cur = eng_buf[(eng_head - 1 + ENERGY_HIST) % ENERGY_HIST];
    double pct = (cur - 1.0) * 100.0;
    char devbuf[16];
    // Format as e.g. "+0.02%" or "-1.34%"
    snprintf(devbuf, sizeof(devbuf), "%+.2f PCT", pct);
    ds(rend, "DRIFT", px, ty + gh + 2, 2, 0xAA, 0xBB, 0xCC);
    Uint8 dcr = (Uint8)(fabs(pct) > 1.0 ? 220 : fabs(pct) > 0.1 ? 180 : 100);
    Uint8 dcg = (Uint8)(fabs(pct) > 1.0 ? 80 : fabs(pct) > 0.1 ? 180 : 220);
    ds(rend, devbuf, px + 40, ty + gh + 2, 2, dcr, dcg, 0x44);
  }
#undef TRACE_Y
  ty += gh + 18;

  // ── KE distribution histogram + Maxwell-Boltzmann curve ───────────────────
  ds(rend, "KE DISTRIBUTION vs MB", px, ty, 2, 0xAA, 0xBB, 0xCC);
  ty += 13;

  int hx = cfg.sim_w + 4, hw = PANEL_W - 8, hh = 85, hy = ty;
  SDL_SetRenderDrawColor(rend, 0x06, 0x08, 0x10, 0xFF);
  SDL_Rect hb = {hx, hy, hw, hh};
  SDL_RenderFillRect(rend, &hb);

  // Grid lines
  SDL_SetRenderDrawBlendMode(rend, SDL_BLENDMODE_BLEND);
  SDL_SetRenderDrawColor(rend, 0x22, 0x33, 0x55, 80);
  for (int f = 1; f <= 3; f++)
    SDL_RenderDrawLine(rend, hx, hy + hh - hh * f / 4, hx + hw,
                       hy + hh - hh * f / 4);
  SDL_SetRenderDrawBlendMode(rend, SDL_BLENDMODE_NONE);

  double bin_w = (double)hw / HIST_BINS;
  for (int b = 0; b < HIST_BINS; b++) {
    int bx = (int)(hx + b * bin_w);
    int bw = (int)bin_w - 1;
    if (bw < 1)
      bw = 1;

    // Measured histogram — cyan bars
    int bar_h = (int)(hist_smooth[b] * (hh - 2));
    if (bar_h > 0) {
      SDL_SetRenderDrawBlendMode(rend, SDL_BLENDMODE_BLEND);
      SDL_SetRenderDrawColor(rend, 0x22, 0xAA, 0xCC, 180);
      SDL_Rect br = {bx, hy + hh - bar_h, bw, bar_h};
      SDL_RenderFillRect(rend, &br);
      SDL_SetRenderDrawBlendMode(rend, SDL_BLENDMODE_NONE);
    }

    // MB theoretical curve — bright yellow line connecting bin tops
    int mb_y = (int)(hy + hh - hist_mb[b] * (hh - 2));
    if (b > 0) {
      int prev_mb_y = (int)(hy + hh - hist_mb[b - 1] * (hh - 2));
      int bx_prev = (int)(hx + (b - 0.5) * bin_w);
      int bx_cur = (int)(hx + (b + 0.5) * bin_w);
      SDL_SetRenderDrawColor(rend, 0xFF, 0xDD, 0x33, 0xFF);
      SDL_RenderDrawLine(rend, bx_prev, prev_mb_y, bx_cur, mb_y);
    }
    (void)mb_y;
  }
  // Draw MB curve endpoint
  {
    int prev_y = (int)(hy + hh - hist_mb[HIST_BINS - 2] * (hh - 2));
    int last_y = (int)(hy + hh - hist_mb[HIST_BINS - 1] * (hh - 2));
    SDL_SetRenderDrawColor(rend, 0xFF, 0xDD, 0x33, 0xFF);
    SDL_RenderDrawLine(rend, (int)(hx + (HIST_BINS - 1.5) * bin_w), prev_y,
                       (int)(hx + (HIST_BINS - 0.5) * bin_w), last_y);
  }

  SDL_SetRenderDrawColor(rend, 0x2A, 0x3A, 0x55, 0xFF);
  SDL_RenderDrawRect(rend, &hb);
  ty += hh + 6;

  // Legend
  SDL_SetRenderDrawColor(rend, 0x22, 0xAA, 0xCC, 0xFF);
  SDL_Rect lc = {px, ty + 2, 8, 8};
  SDL_RenderFillRect(rend, &lc);
  ds(rend, "MEASURED", px + 11, ty, 2, 0x44, 0xCC, 0xEE);
  ty += 14;
  SDL_SetRenderDrawColor(rend, 0xFF, 0xDD, 0x33, 0xFF);
  SDL_RenderDrawLine(rend, px, ty + 4, px + 8, ty + 4);
  ds(rend, "MAXWELL-BOLTZMANN", px + 11, ty, 2, 0xDD, 0xBB, 0x44);
  ty += 20;

  // Species legend
  SDL_SetRenderDrawColor(rend, 0x33, 0xCC, 0xFF, 0xFF);
  SDL_Rect sl = {px, ty + 2, 8, 8};
  SDL_RenderFillRect(rend, &sl);
  ds(rend, "LIGHT", px + 11, ty, 2, 0x33, 0xCC, 0xFF);
  SDL_SetRenderDrawColor(rend, 0xFF, 0x88, 0x33, 0xFF);
  SDL_Rect sh = {px + 60, ty + 2, 8, 8};
  SDL_RenderFillRect(rend, &sh);
  ds(rend, "HEAVY", px + 71, ty, 2, 0xFF, 0x88, 0x33);
  ty += 14;

  // Divider
  SDL_SetRenderDrawColor(rend, 0x2A, 0x3A, 0x55, 0xFF);
  SDL_RenderDrawLine(rend, cfg.sim_w + 4, ty, cfg.sim_w + PANEL_W - 4, ty);
  ty += 6;

  // Sliders
  g_slider_panel_top = ty;
  draw_sliders(rend, ty);

  // Key hints at very bottom
  ds(rend, "R-RAND  M-MONO  V-REC  Q-QUIT", px, cfg.sim_h - 14, 1, 0x55, 0x66,
     0x88);
}

// ── Main
// ──────────────────────────────────────────────────────────────────────
int main(int argc, char **argv) {
  // Parse --config <file> or --write-config <file>
  const char *config_path = "ballsim.cfg";
  for (int i = 1; i < argc; i++) {
    if (!strcmp(argv[i], "--config") && i + 1 < argc) {
      config_path = argv[++i];
    } else if (!strcmp(argv[i], "--write-config")) {
      const char *out = (i + 1 < argc) ? argv[++i] : "ballsim.cfg";
      config_derive();
      write_default_config(out);
      return 0;
    } else if (!strcmp(argv[i], "--help")) {
      printf("Usage: ballsim [--config FILE] [--write-config [FILE]]\n"
             "  --config FILE        load config from FILE (default: "
             "ballsim.cfg)\n"
             "  --write-config FILE  write default config to FILE and exit\n");
      return 0;
    }
  }
  load_config(config_path);
  config_derive();

  srand((unsigned)time(NULL));

  if (SDL_Init(SDL_INIT_VIDEO) != 0) {
    fprintf(stderr, "SDL: %s\n", SDL_GetError());
    return 1;
  }
  SDL_Window *win = SDL_CreateWindow(
      "Parabolic Vessel — Hard Sphere Gas", SDL_WINDOWPOS_CENTERED,
      SDL_WINDOWPOS_CENTERED, CFG_WIN_W, cfg.sim_h, SDL_WINDOW_SHOWN);
  SDL_Renderer *rend = SDL_CreateRenderer(
      win, -1, SDL_RENDERER_ACCELERATED | SDL_RENDERER_PRESENTVSYNC);
  if (!win || !rend) {
    fprintf(stderr, "SDL: %s\n", SDL_GetError());
    return 1;
  }

  int n = cfg.n_particles;
  init_particles(n, true);
  sliders_init();

  Uint32 t_prev = SDL_GetTicks(), t_last_samp = t_prev;
  double accumulator = 0.0;
  const double FIXED_DT = 1.0 / 120.0;
  bool running = true;

  while (running) {
    SDL_Event ev;
    while (SDL_PollEvent(&ev)) {
      if (ev.type == SDL_QUIT)
        running = false;

      if (ev.type == SDL_MOUSEBUTTONDOWN &&
          ev.button.button == SDL_BUTTON_LEFT) {
        int hit = slider_hit(ev.button.x, ev.button.y, g_slider_panel_top);
        if (hit >= 0) {
          active_slider = hit;
          int needs_reset = slider_set(hit, ev.button.x, g_slider_panel_top);
          if (needs_reset) {
            n = (int)cfg.n_particles;
            init_particles(n, true);
          }
        }
      }
      if (ev.type == SDL_MOUSEBUTTONUP && ev.button.button == SDL_BUTTON_LEFT) {
        active_slider = -1;
      }
      if (ev.type == SDL_MOUSEMOTION && active_slider >= 0) {
        int needs_reset =
            slider_set(active_slider, ev.motion.x, g_slider_panel_top);
        if (needs_reset) {
          n = (int)cfg.n_particles;
          init_particles(n, true);
        }
      }

      if (ev.type == SDL_KEYDOWN)
        switch (ev.key.keysym.sym) {
        case SDLK_ESCAPE:
        case SDLK_q:
          running = false;
          break;
        case SDLK_r:
          init_particles(n, false);
          break;
        case SDLK_m:
          init_particles(n, true);
          break;
        case SDLK_EQUALS:
        case SDLK_PLUS:
          if (n < MAX_PARTICLES) {
            n++;
            cfg.n_particles = n;
          }
          break;
        case SDLK_MINUS:
          if (n > 2) {
            n--;
            cfg.n_particles = n;
          }
          break;
        case SDLK_v:
          if (!rec_active) {
            record_start(CFG_WIN_W, cfg.sim_h);
            rec_blink = SDL_GetTicks();
          } else {
            record_stop();
          }
          break;
        }
    }

    Uint32 t_now = SDL_GetTicks();
    double frame_time = (t_now - t_prev) / 1000.0;
    if (frame_time > 0.05)
      frame_time = 0.05; // clamp to avoid spiral of death
    t_prev = t_now;
    accumulator += frame_time;

    // Consume accumulated time in fixed steps
    while (accumulator >= FIXED_DT) {
      step_sim(n, FIXED_DT);
      accumulator -= FIXED_DT;
    }

    if (t_now - t_last_samp >= 80) {
      double ke, pe;
      push_energy(compute_energy(n, &ke, &pe));
      update_histogram(n);
      t_last_samp = t_now;
    }

    SDL_SetRenderDrawColor(rend, 0x0D, 0x0F, 0x1A, 0xFF);
    SDL_RenderClear(rend);
    for (int i = 0; i < n; i++)
      draw_particle(rend, p[i].x, p[i].y, p[i].radius, p[i].r, p[i].g, p[i].b);
    draw_vessel(rend);
    double ke, pe;
    double te = compute_energy(n, &ke, &pe);
    draw_panel(rend, n, te, ke, pe);
    SDL_RenderPresent(rend);
    capture_frame(rend, CFG_WIN_W, cfg.sim_h);
  }

  record_stop();   // flush and close pipe if recording was active

  SDL_DestroyRenderer(rend);
  SDL_DestroyWindow(win);
  SDL_Quit();
  return 0;
}
