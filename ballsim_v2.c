// ballsim_v2.c — Parabolic Vessel · Hard Sphere Gas  (v2.1)
//
// Physics engine: fixed-dt Verlet integration + cell-list O(N) collision
// detection.  Up to 5 000 particles.  OpenMP parallelism over a 9-colour
// checkerboard of cells so collision resolution itself runs multi-core.
//
// Per outer fixed-dt step (120 Hz):
//   for each sub-step (g_substeps per outer dt):
//     1. Verlet integrate under gravity              [O(N), parallel]
//     2. Boundary conditions (walls/ceiling/floor)   [O(N), parallel]
//        — accumulates wall impulse (pressure) and piston work
//     3. Build cell list + colour buckets            [O(N), serial]
//     4. Resolve particle-particle overlaps          [O(N), parallel by colour]
//        repeated g_passes times for dense packings
//
// Analysis column (right of the stats panel):
//   equipartition time series, H-function, barometric density profile,
//   wall pressure + compressibility Z + P-V loop + piston work/cycle,
//   tracer trail + ensemble MSD, mean free path vs Enskog theory.
//
// Modes:
//   ./ballsim_v2                          interactive
//   ./ballsim_v2 --bench N STEPS          headless benchmark (no window)
//   ./ballsim_v2 --screenshot out.bmp K   run K steps, save one frame, exit
//   flags: --config F  --write-config [F]  --substeps K  --passes P
//          --threads T  --seed S  --help
//
// Build (see Makefile):
//   gcc -O3 -std=c11 -fopenmp -o ballsim_v2 ballsim_v2.c $(pkg-config --cflags --libs sdl2) -lm


//
// v1 (event-driven, exact CCD, 500 particles) is kept intact in
//   verlet_collisions_two_gasses_video.c

#include <SDL2/SDL.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#ifdef _OPENMP
#  include <omp.h>
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// ── Compile-time caps ─────────────────────────────────────────────────────────
#define MAX_PARTICLES   5000
#define SERIES_LEN       280   // samples in every time-series ring buffer
#define HIST_BINS         30   // KE-distribution histogram bins
#define PANEL_W          280   // stats panel width
#define ANALYSIS_W       300   // analysis panel width
#define CELL_GRID_MAX    512   // max cells per axis
#define NB_Y              24   // density-profile bins
#define NB_H              40   // H-function bins
#define TRAIL_LEN        240   // tracer trail points
#define SPRITE_S          64   // particle sprite texture size (px)
#define SPRITE_R        30.0   // sphere radius inside the sprite (px)
#define SAMPLE_MS         80   // diagnostics sampling period (wall-clock ms)

// ── Runtime-tunable solver parameters (CLI: --substeps, --passes) ─────────────
static int g_substeps = 8;   // Verlet sub-steps per outer 1/120 s step
static int g_passes   = 4;   // overlap-resolution passes per sub-step

// ── Runtime config ────────────────────────────────────────────────────────────
typedef struct {
    int    n_particles;
    double gravity;
    double radius;
    double restitution;
    double mono_speed;
    double para_a;
    int    sim_w;
    int    sim_h;
    double ceil_amp;
    double ceil_freq;
    double heavy_frac;
    double mass_ratio;
} Config;

static Config cfg = {
    .n_particles = 300,
    .gravity     = 500.0,
    .radius      = 4.0,
    .restitution = 1.0,
    .mono_speed  = 220.0,
    .para_a      = 0.0030,
    .sim_w       = 700,
    .sim_h       = 960,
    .ceil_amp    = 0.0,
    .ceil_freq   = 0.5,
    .heavy_frac  = 0.5,
    .mass_ratio  = 4.0,
};

static double CFG_PARA_BASE;   // y of parabola vertex (floor bottom)
static int    CFG_WIN_W;
static double CFG_ARC;         // floor arc length below the ceiling top
static double CFG_AREA0;       // vessel area with ceiling at its top position
static double CFG_RIM_Y;       // y where the parabola meets the side walls

// ── Parabola vessel ───────────────────────────────────────────────────────────
static inline double para_y(double x) {
    double d = x - cfg.sim_w * 0.5;
    return CFG_PARA_BASE - cfg.para_a * d * d;
}
static inline double para_dydx(double x) {
    return -2.0 * cfg.para_a * (x - cfg.sim_w * 0.5);
}
// Horizontal width of the vessel at height y (for the barometric profile)
static double vessel_width(double y) {
    if (cfg.para_a < 1e-9) return cfg.sim_w;
    double d = (CFG_PARA_BASE - y) / cfg.para_a;
    if (d <= 0) return 0.0;
    double w = 2.0 * sqrt(d);
    return w < cfg.sim_w ? w : cfg.sim_w;
}

static void config_derive(void) {
    CFG_PARA_BASE = (cfg.para_a > 1e-9) ? cfg.sim_h - 80.0 : cfg.sim_h;
    CFG_WIN_W     = cfg.sim_w + PANEL_W + ANALYSIS_W;
    // Floor geometry (clipped at the ceiling's topmost position)
    const int S = 400;
    double top = cfg.radius, arc = 0, area = 0;
    for (int k = 0; k < S; k++) {
        double x0 = (double)k * cfg.sim_w / S, x1 = (double)(k + 1) * cfg.sim_w / S;
        double y0 = para_y(x0), y1 = para_y(x1);
        if (y0 > top && y1 > top) arc += hypot(x1 - x0, y1 - y0);
        double c0 = y0 > top ? y0 : top, c1 = y1 > top ? y1 : top;
        area += 0.5 * ((c0 - top) + (c1 - top)) * (x1 - x0);
    }
    CFG_ARC   = arc;
    CFG_AREA0 = area;
    CFG_RIM_Y = (cfg.para_a > 1e-9) ? para_y(0.0) : CFG_PARA_BASE;
    if (CFG_RIM_Y < top) CFG_RIM_Y = top;
}

static void clamp_config(void) {
    if (cfg.n_particles < 1)             cfg.n_particles = 1;
    if (cfg.n_particles > MAX_PARTICLES) cfg.n_particles = MAX_PARTICLES;
    if (cfg.gravity < 0)      cfg.gravity = 0;
    if (cfg.radius  < 1)      cfg.radius  = 1;
    if (cfg.radius  > 20)     cfg.radius  = 20;
    if (cfg.restitution < 0)  cfg.restitution = 0;
    if (cfg.restitution > 1)  cfg.restitution = 1;
    if (cfg.mono_speed  < 1)  cfg.mono_speed  = 1;
    if (cfg.para_a < 0)       cfg.para_a = 0;
    if (cfg.para_a > 0.01)    cfg.para_a = 0.01;
    if (cfg.sim_w  < 200)     cfg.sim_w  = 200;
    if (cfg.sim_w  > 2560)    cfg.sim_w  = 2560;
    if (cfg.sim_h  < 200)     cfg.sim_h  = 200;
    if (cfg.sim_h  > 1440)    cfg.sim_h  = 1440;
    if (cfg.ceil_amp  < 0)    cfg.ceil_amp  = 0;
    if (cfg.ceil_amp  > 200)  cfg.ceil_amp  = 200;
    if (cfg.ceil_freq < 0)    cfg.ceil_freq = 0;
    if (cfg.heavy_frac < 0)   cfg.heavy_frac = 0;
    if (cfg.heavy_frac > 1)   cfg.heavy_frac = 1;
    if (cfg.mass_ratio < 1)   cfg.mass_ratio = 1;
    if (cfg.mass_ratio > 20)  cfg.mass_ratio = 20;
}

static void load_config(const char *path) {
    FILE *f = fopen(path, "r");
    if (!f) { fprintf(stderr, "Config: cannot open '%s', using defaults.\n", path); return; }
    char line[256];
    while (fgets(line, sizeof(line), f)) {
        char *p = line;
        while (*p == ' ' || *p == '\t') p++;
        if (*p == '#' || *p == '\n' || *p == '\0') continue;
        char key[64]; double val;
        if (sscanf(p, "%63s = %lf", key, &val) != 2) continue;
        if      (!strcmp(key,"n_particles")) cfg.n_particles = (int)val;
        else if (!strcmp(key,"gravity"))     cfg.gravity     = val;
        else if (!strcmp(key,"radius"))      cfg.radius      = val;
        else if (!strcmp(key,"restitution")) cfg.restitution = val;
        else if (!strcmp(key,"mono_speed"))  cfg.mono_speed  = val;
        else if (!strcmp(key,"para_a"))      cfg.para_a      = val;
        else if (!strcmp(key,"sim_w"))       cfg.sim_w       = (int)val;
        else if (!strcmp(key,"sim_h"))       cfg.sim_h       = (int)val;
        else if (!strcmp(key,"ceil_amp"))    cfg.ceil_amp    = val;
        else if (!strcmp(key,"ceil_freq"))   cfg.ceil_freq   = val;
        else if (!strcmp(key,"heavy_frac"))  cfg.heavy_frac  = val;
        else if (!strcmp(key,"mass_ratio"))  cfg.mass_ratio  = val;
        else if (!strcmp(key,"substeps"))    g_substeps      = (int)val;
        else if (!strcmp(key,"passes"))      g_passes        = (int)val;
        else fprintf(stderr, "Config: unknown key '%s'\n", key);
    }
    fclose(f);
    clamp_config();
}

static void write_default_config(const char *path) {
    FILE *f = fopen(path, "w");
    if (!f) { fprintf(stderr, "Config: cannot write '%s'\n", path); return; }
    fprintf(f,
        "# Parabolic Vessel v2.1 — config file\n"
        "sim_w       = %d\n"
        "sim_h       = %d\n"
        "n_particles = %d\n"
        "gravity     = %.1f\n"
        "radius      = %.1f\n"
        "restitution = %.2f\n"
        "mono_speed  = %.1f\n"
        "para_a      = %.4f\n"
        "ceil_amp    = %.1f\n"
        "ceil_freq   = %.2f\n"
        "heavy_frac  = %.2f\n"
        "mass_ratio  = %.1f\n"
        "substeps    = %d        # Verlet sub-steps per 1/120 s\n"
        "passes      = %d        # overlap-resolution passes per sub-step\n",
        cfg.sim_w, cfg.sim_h, cfg.n_particles, cfg.gravity,
        cfg.radius, cfg.restitution, cfg.mono_speed, cfg.para_a,
        cfg.ceil_amp, cfg.ceil_freq, cfg.heavy_frac, cfg.mass_ratio,
        g_substeps, g_passes);
    fclose(f);
    fprintf(stderr, "Config: wrote defaults to '%s'\n", path);
}

// ── Oscillating ceiling ───────────────────────────────────────────────────────
static double sim_time = 0.0;

static inline double ceil_y(double t) {
    if (cfg.ceil_amp < 1e-9) return cfg.radius;
    return cfg.radius + cfg.ceil_amp * (1.0 + sin(2.0 * M_PI * cfg.ceil_freq * t));
}
static inline double ceil_vy(double t) {
    if (cfg.ceil_amp < 1e-9) return 0.0;
    return cfg.ceil_amp * (2.0 * M_PI * cfg.ceil_freq)
           * cos(2.0 * M_PI * cfg.ceil_freq * t);
}
// Instantaneous vessel area and perimeter (geometric, particle radius ignored)
static inline double vessel_area(double t) {
    return CFG_AREA0 - cfg.sim_w * (ceil_y(t) - cfg.radius);
}
static inline double vessel_perim(double t) {
    double side = CFG_RIM_Y - ceil_y(t);
    if (side < 0) side = 0;
    return CFG_ARC + cfg.sim_w + 2.0 * side;
}

// ── Particle ──────────────────────────────────────────────────────────────────
typedef struct {
    double x, y, vx, vy;
    double mass;
    double radius;
    int    species;   // 0 = light, 1 = heavy
    Uint8  r, g, b;   // species colour
    double hue;       // RAINBOW mode
} Particle;

static Particle p[MAX_PARTICLES];
static double   ref_x[MAX_PARTICLES], ref_y[MAX_PARTICLES]; // MSD reference

// ── Time series ring buffer ───────────────────────────────────────────────────
typedef struct { double v[SERIES_LEN]; int head, count; } Series;

static void series_reset(Series *s) { s->head = s->count = 0; }
static void series_push(Series *s, double x) {
    s->v[s->head] = x;
    s->head = (s->head + 1) % SERIES_LEN;
    if (s->count < SERIES_LEN) s->count++;
}
static double series_at(const Series *s, int k) { // k = 0 → newest
    return s->v[(s->head - 1 - k + SERIES_LEN) % SERIES_LEN];
}
static void series_minmax(const Series *s, double *lo, double *hi) {
    if (s->count == 0) { *lo = 0; *hi = 1; return; }
    *lo = *hi = series_at(s, 0);
    for (int k = 1; k < s->count; k++) {
        double v = series_at(s, k);
        if (v < *lo) *lo = v;
        if (v > *hi) *hi = v;
    }
}

// ── Cell list storage (built each sub-step; also read by diagnostics) ───────
static int    g_cell_head[CELL_GRID_MAX * CELL_GRID_MAX];
static int    g_cell_next[MAX_PARTICLES];
static int    g_dirty[MAX_PARTICLES];           // non-empty cells (for clearing)
static int    g_ndirty = 0;
static int    g_bucket[9][MAX_PARTICLES];       // non-empty cells per colour
static int    g_bcount[9];
static int    g_ncx, g_ncy;
static double g_inv_cs;

// ── Diagnostics state ─────────────────────────────────────────────────────────
static Series s_energy;            // E/E0
static double eng_e0 = 0.0;
static double hist_smooth[HIST_BINS], hist_mb[HIST_BINS], hist_ke_max = 1.0;
static long   collision_count = 0;
static double species_ke[2] = {0, 0};
static int    species_n[2]  = {0, 0};

static Series s_eq[2];             // light / heavy KE/N
static Series s_h;                 // H-function (KE marginal)
static double h_range = 1.0;       // KE range for H bins (fixed at reset)
static double h_eq = 0.0;          // equilibrium H for current kT
static double dens_meas[2][NB_Y];  // smoothed particle counts per y-bin

static double imp_acc = 0.0;       // wall impulse since last sample
static double work_acc = 0.0;      // piston work since last sample
static double piston_work_cycle = 0.0, piston_work_last = 0.0;
static double piston_phase_prev = 0.0;
static long   piston_cycles = 0;
static double P_smooth = 0.0, P_ideal = 0.0, Z_meas = 0.0, Z_hend = 0.0;
static double packing = 0.0;
static Series s_P, s_V;            // pressure, volume (area) for P-V loop

static Series s_msd;               // ensemble MSD of heavy particles
static int    tracer_i = 0;
static SDL_FPoint trail[TRAIL_LEN];
static int    trail_head = 0, trail_count = 0;

static double coll_rate = 0.0;     // collisions per particle per second
static double mfp_meas = 0.0, mfp_theory = 0.0, mean_speed = 0.0;

static double samp_prev_time = 0.0;
static long   samp_prev_coll = 0;

static double perf_step_ms = 0.0, perf_fps = 0.0;

static double randf(void) { return (double)rand() / ((double)RAND_MAX + 1.0); }

// ── Colour cycling ────────────────────────────────────────────────────────────
#define N_COLOR_MODES 5
static int color_mode = 0;
static const char *color_mode_names[N_COLOR_MODES] = {
    "SPECIES", "SPEED", "KE", "HEIGHT", "RAINBOW"
};

static void hsv_to_rgb(double h, double s, double v, Uint8 *r, Uint8 *g, Uint8 *b) {
    if (s < 1e-9) { *r = *g = *b = (Uint8)(v * 255); return; }
    int    i = (int)(h / 60.0) % 6;
    double f = h / 60.0 - floor(h / 60.0);
    double pp = v * (1.0 - s), q = v * (1.0 - f * s), t = v * (1.0 - (1.0 - f) * s);
    double rv, gv, bv;
    switch (i) {
        case 0: rv = v;  gv = t;  bv = pp; break;
        case 1: rv = q;  gv = v;  bv = pp; break;
        case 2: rv = pp; gv = v;  bv = t;  break;
        case 3: rv = pp; gv = q;  bv = v;  break;
        case 4: rv = t;  gv = pp; bv = v;  break;
        default: rv = v; gv = pp; bv = q;  break;
    }
    *r = (Uint8)(rv * 255); *g = (Uint8)(gv * 255); *b = (Uint8)(bv * 255);
}

static void particle_rgb(int i, double max_val, Uint8 *r, Uint8 *g, Uint8 *b) {
    switch (color_mode) {
    case 1: {
        double speed = sqrt(p[i].vx * p[i].vx + p[i].vy * p[i].vy);
        double t = (max_val > 0) ? speed / max_val : 0.0;
        if (t > 1) t = 1;
        hsv_to_rgb(240.0 * (1.0 - t), 1.0, 1.0, r, g, b);
        break;
    }
    case 2: {
        double ke = 0.5 * p[i].mass * (p[i].vx * p[i].vx + p[i].vy * p[i].vy);
        double t = (max_val > 0) ? ke / max_val : 0.0;
        if (t > 1) t = 1;
        hsv_to_rgb(240.0 * (1.0 - t), 1.0, 1.0, r, g, b);
        break;
    }
    case 3: {
        double t = p[i].y / cfg.sim_h;
        if (t < 0) t = 0;
        if (t > 1) t = 1;
        hsv_to_rgb(t * 300.0, 1.0, 1.0, r, g, b);
        break;
    }
    case 4:
        hsv_to_rgb(p[i].hue, 1.0, 1.0, r, g, b);
        break;
    default:
        *r = p[i].r; *g = p[i].g; *b = p[i].b;
    }
}

// ── Maxwell-Boltzmann curve (2D: f(KE) ∝ exp(-KE/kT)) ────────────────────────
static void compute_mb_curve(double mean_ke) {
    if (mean_ke < 1.0) mean_ke = 1.0;
    double peak = 0.0, raw[HIST_BINS];
    for (int b = 0; b < HIST_BINS; b++) {
        double ke = (b + 0.5) / HIST_BINS * hist_ke_max;
        raw[b] = exp(-ke / mean_ke);
        if (raw[b] > peak) peak = raw[b];
    }
    for (int b = 0; b < HIST_BINS; b++)
        hist_mb[b] = (peak > 0) ? raw[b] / peak : 0.0;
}

// ── Particle initialisation ───────────────────────────────────────────────────
static void assign_species(int i, bool monoenergetic) {
    int is_heavy  = (randf() < cfg.heavy_frac) ? 1 : 0;
    p[i].species  = is_heavy;
    p[i].mass     = is_heavy ? cfg.mass_ratio : 1.0;
    p[i].radius   = cfg.radius * cbrt(p[i].mass);

    double angle = randf() * 2.0 * M_PI;
    double speed = monoenergetic ? cfg.mono_speed : 80.0 + randf() * 300.0;
    double v     = speed / sqrt(p[i].mass);   // equal-KE start
    p[i].vx = v * cos(angle);
    p[i].vy = v * sin(angle);

    double br = 0.75 + 0.25 * randf();
    if (is_heavy) {
        p[i].r = (Uint8)(255 * br); p[i].g = (Uint8)(110 * br); p[i].b = (Uint8)(20 * br);
    } else {
        p[i].r = (Uint8)(20 * br);  p[i].g = (Uint8)(190 * br); p[i].b = (Uint8)(255 * br);
    }
    p[i].hue = randf() * 360.0;
}

static void msd_reset(int n) {
    for (int i = 0; i < n; i++) { ref_x[i] = p[i].x; ref_y[i] = p[i].y; }
    series_reset(&s_msd);
}

static void reset_diagnostics(int n) {
    series_reset(&s_energy);
    eng_e0 = 0.0;
    for (int b = 0; b < HIST_BINS; b++) hist_smooth[b] = hist_mb[b] = 0.0;
    hist_ke_max = 1.0;
    collision_count = 0;
    sim_time = 0.0;

    series_reset(&s_eq[0]); series_reset(&s_eq[1]);
    series_reset(&s_h);
    for (int s = 0; s < 2; s++)
        for (int b = 0; b < NB_Y; b++) dens_meas[s][b] = 0.0;

    imp_acc = work_acc = 0.0;
    piston_work_cycle = piston_work_last = 0.0;
    piston_phase_prev = 0.0;
    piston_cycles = 0;
    P_smooth = P_ideal = Z_meas = Z_hend = packing = 0.0;
    series_reset(&s_P); series_reset(&s_V);

    msd_reset(n);
    tracer_i = 0;
    for (int i = 0; i < n; i++) if (p[i].species == 1) { tracer_i = i; break; }
    trail_head = trail_count = 0;

    coll_rate = mfp_meas = mfp_theory = mean_speed = 0.0;
    samp_prev_time = 0.0;
    samp_prev_coll = 0;

    // H-function bin range: 8 kT covers 99.97 % of an exponential distribution
    double ke_sum = 0;
    for (int i = 0; i < n; i++)
        ke_sum += 0.5 * p[i].mass * (p[i].vx * p[i].vx + p[i].vy * p[i].vy);
    h_range = 8.0 * (n > 0 ? ke_sum / n : 1.0);
    if (h_range < 1.0) h_range = 1.0;
}

static void init_particles(int n, bool monoenergetic) {
    double max_r = cfg.radius;
    for (int i = 0; i < n; i++) {
        assign_species(i, monoenergetic);
        if (p[i].radius > max_r) max_r = p[i].radius;
    }
    // Grid placement, stacked from the FLOOR upward and sized for the largest
    // particle.  Starting at the bottom keeps the initial potential energy
    // small, so the temperature is set by mono_speed rather than by the drop
    // (a 950 px fall at g=500 would add 475 000 KE per unit mass), and keeps
    // the lattice clear of the piston at its lowest point.
    double gap    = 2.0;
    double cell   = 2.0 * max_r + gap;
    double y_min  = cfg.radius + 2.0 * cfg.ceil_amp + max_r + gap; // below piston
    double y      = CFG_PARA_BASE - max_r - gap;                    // first row
    int placed = 0;
    while (placed < n && y >= y_min) {
        // Usable half-width at this row: centres must satisfy para_y(x) ≥ y+max_r+gap
        double hw;
        if (cfg.para_a > 1e-9) {
            double d = (CFG_PARA_BASE - y - max_r - gap) / cfg.para_a;
            hw = (d > 0) ? sqrt(d) : 0.0;
        } else hw = cfg.sim_w * 0.5;
        double hw_wall = cfg.sim_w * 0.5 - max_r - gap;
        if (hw > hw_wall) hw = hw_wall;
        int cols = (int)(2.0 * hw / cell);
        if (cols >= 1) {
            double x0 = cfg.sim_w * 0.5 - cols * cell * 0.5 + cell * 0.5;
            for (int c = 0; c < cols && placed < n; c++, placed++) {
                p[placed].x = x0 + c * cell + (randf() - 0.5) * gap;
                p[placed].y = y             + (randf() - 0.5) * gap;
            }
        }
        y -= cell;
    }
    // Vessel full: scatter the remainder at random; the solver separates them
    for (; placed < n; placed++) {
        p[placed].x = max_r + gap + randf() * (cfg.sim_w - 2 * (max_r + gap));
        p[placed].y = y_min + randf() * (CFG_PARA_BASE - max_r - gap - y_min);
    }
    for (int i = 0; i < n; i++) {
        double ri = p[i].radius;
        if (p[i].x < ri + 1)             p[i].x = ri + 1;
        if (p[i].x > cfg.sim_w - ri - 1) p[i].x = cfg.sim_w - ri - 1;
        if (p[i].y < ri + 1)             p[i].y = ri + 1;
        if (p[i].y > cfg.sim_h - ri - 1) p[i].y = cfg.sim_h - ri - 1;
    }
    reset_diagnostics(n);
}

// Spawn one extra particle (+ key) just under the ceiling
static void spawn_particle(int i) {
    assign_species(i, true);
    double ri = p[i].radius;
    p[i].x = ri + 2 + randf() * (cfg.sim_w - 2 * ri - 4);
    p[i].y = ceil_y(sim_time) + ri + 2 + randf() * 30.0;
    ref_x[i] = p[i].x; ref_y[i] = p[i].y;
}

// ── Energy ────────────────────────────────────────────────────────────────────
static double compute_energy(int n, double *ke_out, double *pe_out) {
    double ke = 0, pe = 0;
    species_ke[0] = species_ke[1] = 0;
    species_n[0]  = species_n[1]  = 0;
    for (int i = 0; i < n; i++) {
        double pke = 0.5 * p[i].mass * (p[i].vx * p[i].vx + p[i].vy * p[i].vy);
        ke += pke;
        pe += p[i].mass * cfg.gravity * (cfg.sim_h - p[i].radius - p[i].y);
        species_ke[p[i].species] += pke;
        species_n [p[i].species]++;
    }
    if (ke_out) *ke_out = ke;
    if (pe_out) *pe_out = pe;
    return ke + pe;
}

static void push_energy(double e) {
    if (eng_e0 == 0.0) eng_e0 = (e != 0.0) ? e : 1.0;
    series_push(&s_energy, e / eng_e0);
}

static void update_histogram(int n) {
    double max_ke = 1.0;
    for (int i = 0; i < n; i++) {
        double ke = 0.5 * p[i].mass * (p[i].vx * p[i].vx + p[i].vy * p[i].vy);
        if (ke > max_ke) max_ke = ke;
    }
    hist_ke_max += (max_ke * 1.1 - hist_ke_max) * 0.05;
    double raw[HIST_BINS] = {0};
    for (int i = 0; i < n; i++) {
        double ke = 0.5 * p[i].mass * (p[i].vx * p[i].vx + p[i].vy * p[i].vy);
        int b = (int)(ke / hist_ke_max * HIST_BINS);
        if (b >= HIST_BINS) b = HIST_BINS - 1;
        raw[b] += 1.0;
    }
    double peak = 0;
    for (int b = 0; b < HIST_BINS; b++) if (raw[b] > peak) peak = raw[b];
    if (peak > 0) for (int b = 0; b < HIST_BINS; b++) raw[b] /= peak;
    for (int b = 0; b < HIST_BINS; b++)
        hist_smooth[b] = hist_smooth[b] * 0.85 + raw[b] * 0.15;
    double ke_total = 0;
    for (int i = 0; i < n; i++)
        ke_total += 0.5 * p[i].mass * (p[i].vx * p[i].vx + p[i].vy * p[i].vy);
    compute_mb_curve(n > 0 ? ke_total / n : 1.0);
}

// ── H-function on the KE marginal ─────────────────────────────────────────────
// H = Σ_b p_b ln p_b over fixed KE bins; decreases toward the exponential
// (Maxwell-Boltzmann) value as the gas thermalises (Boltzmann H-theorem).
static double compute_H(int n) {
    int cnt[NB_H] = {0};
    for (int i = 0; i < n; i++) {
        double ke = 0.5 * p[i].mass * (p[i].vx * p[i].vx + p[i].vy * p[i].vy);
        int b = (int)(ke / h_range * NB_H);
        if (b >= NB_H) b = NB_H - 1;
        if (b < 0) b = 0;
        cnt[b]++;
    }
    double H = 0;
    for (int b = 0; b < NB_H; b++)
        if (cnt[b] > 0) { double pb = (double)cnt[b] / n; H += pb * log(pb); }
    return H;
}
static double compute_H_eq(double kT) {
    if (kT < 1e-9) return 0.0;
    double H = 0, bw = h_range / NB_H;
    for (int b = 0; b < NB_H; b++) {
        double a0 = b * bw, a1 = (b + 1) * bw;
        double pb = (b < NB_H - 1) ? exp(-a0 / kT) - exp(-a1 / kT) : exp(-a0 / kT);
        if (pb > 0) H += pb * log(pb);
    }
    return H;
}

// ── Diagnostics sample (called every SAMPLE_MS of wall-clock, or every
//    10 steps in headless mode) ─────────────────────────────────────────────────
static void sample_diagnostics(int n) {
    double ke, pe;
    double te = compute_energy(n, &ke, &pe);
    push_energy(te);
    update_histogram(n);
    double kT = (n > 0) ? ke / n : 0.0;

    // Equipartition series
    series_push(&s_eq[0], species_n[0] > 0 ? species_ke[0] / species_n[0] : 0.0);
    series_push(&s_eq[1], species_n[1] > 0 ? species_ke[1] / species_n[1] : 0.0);

    // H-function
    series_push(&s_h, compute_H(n));
    h_eq = compute_H_eq(kT);

    // Density profile (counts per y-bin between ceiling and floor bottom)
    {
        double y_top = ceil_y(sim_time), y_bot = CFG_PARA_BASE;
        double span = y_bot - y_top; if (span < 1) span = 1;
        double cnt[2][NB_Y] = {{0}};
        for (int i = 0; i < n; i++) {
            int b = (int)((p[i].y - y_top) / span * NB_Y);
            if (b < 0) b = 0;
            if (b >= NB_Y) b = NB_Y - 1;
            cnt[p[i].species][b] += 1.0;
        }
        for (int s = 0; s < 2; s++)
            for (int b = 0; b < NB_Y; b++)
                dens_meas[s][b] = dens_meas[s][b] * 0.9 + cnt[s][b] * 0.1;
    }

    // Pressure, compressibility, collision rate
    double dts = sim_time - samp_prev_time;
    double area = vessel_area(sim_time);
    if (dts > 1e-9) {
        double P_inst = imp_acc / (dts * vessel_perim(sim_time));
        P_smooth = (P_smooth == 0.0) ? P_inst : 0.75 * P_smooth + 0.25 * P_inst;
        double dc = (double)(collision_count - samp_prev_coll);
        double rate_inst = (n > 0) ? 2.0 * dc / (dts * n) : 0.0;
        coll_rate = (coll_rate == 0.0) ? rate_inst : 0.75 * coll_rate + 0.25 * rate_inst;
    }
    imp_acc = 0.0;
    work_acc = 0.0;
    samp_prev_time = sim_time;
    samp_prev_coll = collision_count;

    P_ideal = (area > 0) ? n * kT / area : 0.0;           // 2D ideal gas: PA = NkT
    Z_meas  = (P_ideal > 0) ? P_smooth / P_ideal : 0.0;
    double area_p = 0, r_sum = 0;
    for (int i = 0; i < n; i++) { area_p += M_PI * p[i].radius * p[i].radius; r_sum += p[i].radius; }
    packing = (area > 0) ? area_p / area : 0.0;
    double eta = packing;
    Z_hend  = (eta < 1) ? (1.0 + eta * eta / 8.0) / ((1.0 - eta) * (1.0 - eta)) : 0.0; // Henderson 1975
    series_push(&s_P, P_smooth);
    series_push(&s_V, area);

    // Mean speed, mean free path (Enskog-corrected hard discs)
    double vsum = 0, r2_sum = 0;
    for (int i = 0; i < n; i++) {
        vsum   += sqrt(p[i].vx * p[i].vx + p[i].vy * p[i].vy);
        r2_sum += p[i].radius * p[i].radius;
    }
    mean_speed = (n > 0) ? vsum / n : 0.0;
    mfp_meas   = (coll_rate > 0) ? mean_speed / coll_rate : 0.0;
    {
        // Under gravity the gas is a thin atmosphere, so N/A is not the density a
        // particle sees.  Use the number-weighted mean density over the density-
        // profile bins:  n_eff = Σ_b c_b² / (N·A_b), A_b = width(y_b)·Δy.
        // (Bins are ~26 000 px², so excluded-volume bias is negligible.)
        double y_top = ceil_y(sim_time), span = CFG_PARA_BASE - y_top;
        if (span < 1) span = 1;
        double dy_b = span / NB_Y, num = 0, tot = 0;
        for (int b = 0; b < NB_Y; b++) {
            double c  = dens_meas[0][b] + dens_meas[1][b];
            double Ab = vessel_width(y_top + (b + 0.5) * dy_b) * dy_b;
            if (Ab > 1) num += c * c / Ab;
            tot += c;
        }
        double n_eff = (tot > 0) ? num / tot : 0.0;
        double sigma = (n > 0) ? 2.0 * r_sum / n : 1.0;
        double eta_e = (n > 0) ? n_eff * M_PI * r2_sum / n : 0.0;   // local packing
        if (eta_e > 0.9) eta_e = 0.9;
        double chi   = (1.0 - 7.0 * eta_e / 16.0) / ((1.0 - eta_e) * (1.0 - eta_e));
        // 2D hard discs: a centre collides when it passes within σ on EITHER side
        // of a target, so the cross-section is 2σ:  λ = 1/(2√2 n σ χ)
        mfp_theory   = (n_eff > 0) ? 1.0 / (2.0 * sqrt(2.0) * n_eff * sigma * chi) : 0.0;
    }

    // Ensemble MSD of heavy particles (all particles if no heavies)
    {
        double acc = 0; int cnt = 0;
        for (int i = 0; i < n; i++) {
            if (species_n[1] > 0 && p[i].species != 1) continue;
            double dx = p[i].x - ref_x[i], dy = p[i].y - ref_y[i];
            acc += dx * dx + dy * dy; cnt++;
        }
        series_push(&s_msd, cnt > 0 ? acc / cnt : 0.0);
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
// ── Physics: cell-list Verlet MD with 9-colour parallel resolution ───────────
// ═══════════════════════════════════════════════════════════════════════════════
//
// Cell size ≥ 2·max_r, so overlapping pairs always lie in the same or adjacent
// cells.  Cells are coloured by (cx mod 3, cy mod 3): two cells of one colour
// are ≥ 3 apart, so their 3×3 neighbourhoods are disjoint and all cells of a
// colour can be resolved concurrently without locks.


static void cell_list_init(void) {
    memset(g_cell_head, 0xFF, sizeof(g_cell_head));   // all -1
    g_ndirty = 0;
}

static void build_cells(int n) {
    double max_r = cfg.radius * cbrt(cfg.mass_ratio > 1.0 ? cfg.mass_ratio : 1.0);
    double cs    = 2.0 * max_r + 0.5;
    g_inv_cs = 1.0 / cs;
    g_ncx = (int)(cfg.sim_w * g_inv_cs) + 2;
    g_ncy = (int)(cfg.sim_h * g_inv_cs) + 2;
    if (g_ncx > CELL_GRID_MAX) g_ncx = CELL_GRID_MAX;
    if (g_ncy > CELL_GRID_MAX) g_ncy = CELL_GRID_MAX;

    for (int k = 0; k < g_ndirty; k++) g_cell_head[g_dirty[k]] = -1;
    g_ndirty = 0;
    for (int c = 0; c < 9; c++) g_bcount[c] = 0;

    for (int i = 0; i < n; i++) {
        int cx = (int)(p[i].x * g_inv_cs);
        int cy = (int)(p[i].y * g_inv_cs);
        if (cx < 0) cx = 0;
        if (cx >= g_ncx) cx = g_ncx - 1;
        if (cy < 0) cy = 0;
        if (cy >= g_ncy) cy = g_ncy - 1;
        int cidx = cy * g_ncx + cx;
        if (g_cell_head[cidx] == -1) {
            g_dirty[g_ndirty++] = cidx;
            int col = (cx % 3) + 3 * (cy % 3);
            g_bucket[col][g_bcount[col]++] = cidx;
        }
        g_cell_next[i]    = g_cell_head[cidx];
        g_cell_head[cidx] = i;
    }
}

// Energy-consistent vertical move.  Any position correction outside the
// integrator (wall overshoot, overlap push) changes gravitational PE by
// m·g·Δy; without compensation that is a steady energy injection (measured:
// +0.04 %/s at N=1000, g=500).  Rescaling |v| so that v² += 2·g·Δy makes the
// shift energetically identical to free fall over the same height.
// y grows downward, so Δy > 0 is a drop (gains KE).
static inline void move_y(Particle *q, double newy) {
    double dy = newy - q->y;
    q->y = newy;
    if (cfg.gravity <= 0.0 || dy == 0.0) return;
    double v2  = q->vx * q->vx + q->vy * q->vy;
    double v2n = v2 + 2.0 * cfg.gravity * dy;
    if (v2n < 0) v2n = 0;
    if (v2 > 1e-18) { double s = sqrt(v2n / v2); q->vx *= s; q->vy *= s; }
}

// Boundary conditions.  Overshoot past a wall is mirrored back (the particle
// is placed where it would have been had it reflected at contact) and the
// velocity reflected.  Accumulates |Δp| at every wall (pressure) and the
// energy the moving piston transfers to particles (work).
static void apply_boundaries(int n, double t_abs) {
    const double cy  = ceil_y(t_abs);
    const double cvy = ceil_vy(t_abs);
    const double e   = cfg.restitution;
    double imp = 0.0, work = 0.0;

#ifdef _OPENMP
#   pragma omp parallel for schedule(static) reduction(+:imp,work)
#endif
    for (int i = 0; i < n; i++) {
        double ri = p[i].radius, m = p[i].mass;

        // Ceiling / piston
        if (p[i].y - ri < cy) {
            double pen = (cy + ri) - p[i].y;                 // penetration depth
            double rel = p[i].vy - cvy;
            if (rel < 0) {
                double ke0 = 0.5 * m * (p[i].vx * p[i].vx + p[i].vy * p[i].vy);
                double vy1 = -rel * e + cvy;
                imp  += m * fabs(vy1 - p[i].vy);
                p[i].vy = vy1;
                work += 0.5 * m * (p[i].vx * p[i].vx + p[i].vy * p[i].vy) - ke0;
                move_y(&p[i], cy + ri + pen * e + 1e-4);     // mirror
            } else {
                move_y(&p[i], cy + ri + 1e-4);               // already separating: clamp
            }
        }
        // Side walls (above the parabola rim only) — horizontal, no PE change
        if (cfg.para_a < 1e-9 || p[i].y < para_y(ri)) {
            if (p[i].x - ri < 0) {
                double pen = ri - p[i].x;
                if (p[i].vx < 0) { imp += m * (1 + e) * (-p[i].vx); p[i].vx = -p[i].vx * e; p[i].x = ri + pen * e + 1e-4; }
                else p[i].x = ri + 1e-4;
            }
            if (p[i].x + ri > cfg.sim_w) {
                double pen = p[i].x + ri - cfg.sim_w;
                if (p[i].vx > 0) { imp += m * (1 + e) * p[i].vx; p[i].vx = -p[i].vx * e; p[i].x = cfg.sim_w - ri - pen * e - 1e-4; }
                else p[i].x = cfg.sim_w - ri - 1e-4;
            }
        }
        // Floor
        if (cfg.para_a > 1e-9) {
            double floor_y = para_y(p[i].x);
            if (p[i].y + ri > floor_y) {
                double pen = p[i].y + ri - floor_y;
                double dydx = para_dydx(p[i].x);
                double nx = -dydx, ny = 1.0, nl = sqrt(nx * nx + ny * ny);
                nx /= nl; ny /= nl;
                double vn = p[i].vx * nx + p[i].vy * ny;
                if (vn > 0) {
                    imp += m * (1 + e) * vn;
                    p[i].vx = (p[i].vx - 2.0 * vn * nx) * e;
                    p[i].vy = (p[i].vy - 2.0 * vn * ny) * e;
                    // mirror along the (near-vertical) normal: vertical part only,
                    // the horizontal component of a shallow-slope mirror is tiny
                    move_y(&p[i], floor_y - ri - pen * e * ny - 1e-4);
                } else {
                    move_y(&p[i], floor_y - ri - 1e-4);
                }
            }
        } else {
            double floor_y = CFG_PARA_BASE - ri;
            if (p[i].y > floor_y) {
                double pen = p[i].y - floor_y;
                if (p[i].vy > 0) {
                    imp += m * (1 + e) * p[i].vy; p[i].vy = -p[i].vy * e;
                    move_y(&p[i], floor_y - pen * e - 1e-4);
                } else {
                    move_y(&p[i], floor_y - 1e-4);
                }
            }
        }
    }
    imp_acc += imp;
    work_acc += work;
    piston_work_cycle += work;
}

static inline void resolve_pair(int i, int j, long *cc) {
    double dx  = p[j].x - p[i].x, dy = p[j].y - p[i].y;
    double d2  = dx * dx + dy * dy;
    double tch = p[i].radius + p[j].radius;
    if (d2 >= tch * tch || d2 < 1e-20) return;
    double d  = sqrt(d2);
    double nx = dx / d, ny = dy / d;
    double mi = p[i].mass, mj = p[j].mass;

    // Elastic impulse (with restitution) if approaching
    double dvn = (p[j].vx - p[i].vx) * nx + (p[j].vy - p[i].vy) * ny;
    if (dvn < 0) {
        double J = -(1.0 + cfg.restitution) * mi * mj / (mi + mj) * dvn;
        p[i].vx -= (J / mi) * nx;  p[i].vy -= (J / mi) * ny;
        p[j].vx += (J / mj) * nx;  p[j].vy += (J / mj) * ny;
        (*cc)++;
    }
    // Positional separation, split inversely by mass (energy-consistent)
    double overlap = tch - d;
    if (overlap <= 0) return;
    double push = overlap * 0.5 + 1e-5, tot = mi + mj;
    p[i].x -= push * (mj / tot) * nx;  move_y(&p[i], p[i].y - push * (mj / tot) * ny);
    p[j].x += push * (mi / tot) * nx;  move_y(&p[j], p[j].y + push * (mi / tot) * ny);
}

static void resolve_collisions(void) {
    long cc = 0;
    for (int col = 0; col < 9; col++) {
        const int  nc    = g_bcount[col];
        const int *cells = g_bucket[col];
#ifdef _OPENMP
#       pragma omp parallel for schedule(dynamic, 16) reduction(+:cc)
#endif
        for (int k = 0; k < nc; k++) {
            int cidx = cells[k];
            int cx = cidx % g_ncx, cy = cidx / g_ncx;
            for (int i = g_cell_head[cidx]; i >= 0; i = g_cell_next[i]) {
                for (int dcy = -1; dcy <= 1; dcy++) {
                    int ny2 = cy + dcy;
                    if (ny2 < 0 || ny2 >= g_ncy) continue;
                    for (int dcx = -1; dcx <= 1; dcx++) {
                        int nx2 = cx + dcx;
                        if (nx2 < 0 || nx2 >= g_ncx) continue;
                        int nidx = ny2 * g_ncx + nx2;
                        for (int j = g_cell_head[nidx]; j >= 0; j = g_cell_next[j]) {
                            if (j <= i) continue;        // each pair once
                            resolve_pair(i, j, &cc);
                        }
                    }
                }
            }
        }
    }
    collision_count += cc;
}

static void step_sim(int n, double dt) {
    double sub_dt = dt / g_substeps;
    const double g = cfg.gravity;

    for (int sub = 0; sub < g_substeps; sub++) {
        double t_end = sim_time + (sub + 1) * sub_dt;

#ifdef _OPENMP
#       pragma omp parallel for schedule(static)
#endif
        for (int i = 0; i < n; i++) {
            p[i].x  += p[i].vx * sub_dt;
            p[i].y  += p[i].vy * sub_dt + 0.5 * g * sub_dt * sub_dt;
            p[i].vy += g * sub_dt;
        }

        apply_boundaries(n, t_end);
        build_cells(n);

        const double cy2 = ceil_y(t_end);
        for (int pass = 0; pass < g_passes; pass++) {
            resolve_collisions();
            // Re-clamp to hard walls: position corrections can push out
#ifdef _OPENMP
#           pragma omp parallel for schedule(static)
#endif
            for (int i = 0; i < n; i++) {
                double ri = p[i].radius;
                if (p[i].x < ri)             p[i].x = ri;
                if (p[i].x > cfg.sim_w - ri) p[i].x = cfg.sim_w - ri;
                if (p[i].y - ri < cy2)       move_y(&p[i], cy2 + ri);
                if (cfg.para_a < 1e-9) {
                    if (p[i].y > CFG_PARA_BASE - ri) move_y(&p[i], CFG_PARA_BASE - ri);
                } else {
                    double fy = para_y(p[i].x);
                    if (p[i].y + ri > fy) move_y(&p[i], fy - ri);
                }
            }
        }
    }
    sim_time += dt;

    // Piston cycle bookkeeping (work per full period)
    if (cfg.ceil_amp > 1e-9 && cfg.ceil_freq > 0) {
        double phase = fmod(cfg.ceil_freq * sim_time, 1.0);
        if (phase < piston_phase_prev) {
            piston_work_last  = piston_work_cycle;
            piston_work_cycle = 0.0;
            piston_cycles++;
        }
        piston_phase_prev = phase;
    }
}

// Diagnostic: largest residual overlap as a fraction of touch distance
static double max_overlap_fraction(int n) {
    build_cells(n);
    double worst = 0;
    for (int c = 0; c < g_ndirty; c++) {
        int cidx = g_dirty[c];
        int cx = cidx % g_ncx, cy = cidx / g_ncx;
        for (int i = g_cell_head[cidx]; i >= 0; i = g_cell_next[i])
            for (int dcy = -1; dcy <= 1; dcy++) for (int dcx = -1; dcx <= 1; dcx++) {
                int nx2 = cx + dcx, ny2 = cy + dcy;
                if (nx2 < 0 || nx2 >= g_ncx || ny2 < 0 || ny2 >= g_ncy) continue;
                for (int j = g_cell_head[ny2 * g_ncx + nx2]; j >= 0; j = g_cell_next[j]) {
                    if (j <= i) continue;
                    double dx = p[j].x - p[i].x, dy = p[j].y - p[i].y;
                    double tch = p[i].radius + p[j].radius;
                    double f = (tch - sqrt(dx * dx + dy * dy)) / tch;
                    if (f > worst) worst = f;
                }
            }
    }
    return worst;
}

// ── Headless benchmark ────────────────────────────────────────────────────────
static double now_s(void) {
    struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

static int run_bench(int n, int steps) {
    const double FIXED_DT = 1.0 / 120.0;
    cfg.n_particles = n;
    clamp_config(); config_derive();
    n = cfg.n_particles;
    cell_list_init();
    init_particles(n, true);
    double ke, pe, e0 = compute_energy(n, &ke, &pe);
    int threads = 1;
#ifdef _OPENMP
    threads = omp_get_max_threads();
#endif
    printf("bench: N=%d steps=%d substeps=%d passes=%d threads=%d radius=%.1f\n",
           n, steps, g_substeps, g_passes, threads, cfg.radius);
    double t0 = now_s();
    for (int s = 0; s < steps; s++) {
        step_sim(n, FIXED_DT);
        if (s % 10 == 9) sample_diagnostics(n);
    }
    double wall = now_s() - t0;
    double e1 = compute_energy(n, &ke, &pe);
    double sim_s = steps * FIXED_DT;
    printf("  wall %.2f s | %.0f steps/s | %.2f ms/step | realtime x%.2f\n",
           wall, steps / wall, 1000.0 * wall / steps, sim_s / wall);
    printf("  energy drift %+.4f %% | collisions %ld (%.1f per particle per s)\n",
           100.0 * (e1 / e0 - 1.0), collision_count,
           2.0 * collision_count / (n * sim_s));
    printf("  packing %.3f | Z meas %.3f | Z Henderson %.3f | mfp meas %.1f / theory %.1f px\n",
           packing, Z_meas, Z_hend, mfp_meas, mfp_theory);
    printf("  max residual overlap %.2f %% of touch distance\n",
           100.0 * max_overlap_fraction(n));
    return 0;
}

// ═══════════════════════════════════════════════════════════════════════════════
// ── Video recording (ffmpeg pipe) ─────────────────────────────────────────────
// ═══════════════════════════════════════════════════════════════════════════════
#define REC_FPS 60
static FILE  *rec_pipe    = NULL;
static bool   rec_active  = false;
static Uint32 rec_blink   = 0;
static long   rec_frames  = 0;
static Uint8 *rec_pixels  = NULL;
static int    rec_w = 0, rec_h = 0;   // renderer output size (HiDPI-aware)
static char   rec_filename[64];
static Uint32 rec_wall_next = 0;

static bool record_start(SDL_Renderer *rend) {
    if (rec_active) return true;
    SDL_GetRendererOutputSize(rend, &rec_w, &rec_h);
    rec_w &= ~1; rec_h &= ~1;          // yuv420p needs even dimensions
    time_t now = time(NULL);
    strftime(rec_filename, sizeof(rec_filename),
             "ballsim_v2_%Y%m%d_%H%M%S.mp4", localtime(&now));
    char cmd[512];
    snprintf(cmd, sizeof(cmd),
        "ffmpeg -y -f rawvideo -pixel_format rgb24"
        " -video_size %dx%d -framerate %d -i pipe:0"
        " -c:v libx264 -preset fast -crf 18 -pix_fmt yuv420p"
        " \"%s\" 2>/dev/null",
        rec_w, rec_h, REC_FPS, rec_filename);
    rec_pipe = popen(cmd, "w");
    if (!rec_pipe) { fprintf(stderr, "Record: ffmpeg not found?\n"); return false; }
    free(rec_pixels);
    rec_pixels = (Uint8 *)malloc((size_t)rec_w * rec_h * 3);
    if (!rec_pixels) { pclose(rec_pipe); rec_pipe = NULL; return false; }
    rec_active = true; rec_frames = 0; rec_wall_next = SDL_GetTicks();
    fprintf(stderr, "Record: started %dx%d → %s\n", rec_w, rec_h, rec_filename);
    return true;
}
static void record_stop(void) {
    if (!rec_active) return;
    rec_active = false;
    if (rec_pipe) { pclose(rec_pipe); rec_pipe = NULL; }
    fprintf(stderr, "Record: %ld frames → %s\n", rec_frames, rec_filename);
}
static void capture_frame(SDL_Renderer *rend) {
    if (!rec_active || !rec_pipe || !rec_pixels) return;
    Uint32 now = SDL_GetTicks();
    if (now < rec_wall_next) return;
    SDL_Rect r = {0, 0, rec_w, rec_h};
    SDL_RenderReadPixels(rend, &r, SDL_PIXELFORMAT_RGB24, rec_pixels, rec_w * 3);
    fwrite(rec_pixels, 1, (size_t)rec_w * rec_h * 3, rec_pipe);
    rec_frames++;
    rec_wall_next += (Uint32)(1000 / REC_FPS);
    if (now > rec_wall_next + (Uint32)(1000 / REC_FPS)) rec_wall_next = now;
}

// ═══════════════════════════════════════════════════════════════════════════════
// ── Drawing ───────────────────────────────────────────────────────────────────
// ═══════════════════════════════════════════════════════════════════════════════

// Particle sprite: white shaded disc with anti-aliased edge; tinted per
// particle through vertex colour, all particles drawn in ONE RenderGeometry
// call (fallback: RenderCopy per particle on SDL < 2.0.18).
static SDL_Texture *g_sprite = NULL;

static SDL_Texture *make_sprite(SDL_Renderer *rend) {
    static Uint8 px[SPRITE_S * SPRITE_S * 4];
    double c = (SPRITE_S - 1) * 0.5;
    for (int y = 0; y < SPRITE_S; y++)
        for (int x = 0; x < SPRITE_S; x++) {
            double dx = x - c, dy = y - c, d = sqrt(dx * dx + dy * dy);
            double a = (SPRITE_R + 0.75 - d) / 1.5;
            if (a < 0) a = 0;
            if (a > 1) a = 1;
            double u = d / SPRITE_R;
            double br = 1.0 - 0.5 * u * u;
            double sx = dx + 0.35 * SPRITE_R, sy = dy + 0.35 * SPRITE_R;
            br += 0.35 * exp(-(sx * sx + sy * sy) / (0.18 * SPRITE_R * 0.18 * SPRITE_R * 4));
            if (br > 1) br = 1;
            if (br < 0) br = 0;
            Uint8 *q = &px[(y * SPRITE_S + x) * 4];
            q[0] = q[1] = q[2] = (Uint8)(br * 255);
            q[3] = (Uint8)(a * 255);
        }
    SDL_Texture *t = SDL_CreateTexture(rend, SDL_PIXELFORMAT_RGBA32,
                                       SDL_TEXTUREACCESS_STATIC, SPRITE_S, SPRITE_S);
    if (!t) return NULL;
    SDL_UpdateTexture(t, NULL, px, SPRITE_S * 4);
    SDL_SetTextureBlendMode(t, SDL_BLENDMODE_BLEND);
    return t;
}

static void draw_particles(SDL_Renderer *rend, int n, double max_val) {
    const float k = (float)((SPRITE_S * 0.5) / SPRITE_R);   // quad half-size / radius
#if SDL_VERSION_ATLEAST(2, 0, 18)
    static SDL_Vertex V[4 * MAX_PARTICLES];
    static int        I[6 * MAX_PARTICLES];
    for (int i = 0; i < n; i++) {
        Uint8 r, g, b;
        particle_rgb(i, max_val, &r, &g, &b);
        float h = (float)p[i].radius * k;
        float x = (float)p[i].x, y = (float)p[i].y;
        SDL_Color col = {r, g, b, 255};
        SDL_Vertex *v = &V[4 * i];
        v[0] = (SDL_Vertex){{x - h, y - h}, col, {0, 0}};
        v[1] = (SDL_Vertex){{x + h, y - h}, col, {1, 0}};
        v[2] = (SDL_Vertex){{x + h, y + h}, col, {1, 1}};
        v[3] = (SDL_Vertex){{x - h, y + h}, col, {0, 1}};
        int *ix = &I[6 * i];
        ix[0] = 4*i; ix[1] = 4*i+1; ix[2] = 4*i+2;
        ix[3] = 4*i; ix[4] = 4*i+2; ix[5] = 4*i+3;
    }
    if (n > 0) SDL_RenderGeometry(rend, g_sprite, V, 4 * n, I, 6 * n);
#else
    for (int i = 0; i < n; i++) {
        Uint8 r, g, b;
        particle_rgb(i, max_val, &r, &g, &b);
        SDL_SetTextureColorMod(g_sprite, r, g, b);
        float h = (float)p[i].radius * k;
        SDL_Rect dst = {(int)(p[i].x - h), (int)(p[i].y - h), (int)(2 * h), (int)(2 * h)};
        SDL_RenderCopy(rend, g_sprite, NULL, &dst);
    }
#endif
}

static void draw_tracer(SDL_Renderer *rend, int n) {
    if (tracer_i >= n) return;
    SDL_SetRenderDrawBlendMode(rend, SDL_BLENDMODE_BLEND);
    // Trail, fading from old to new
    for (int k = trail_count - 1; k >= 1; k--) {
        SDL_FPoint a = trail[(trail_head - 1 - k + TRAIL_LEN) % TRAIL_LEN];
        SDL_FPoint b = trail[(trail_head - k + TRAIL_LEN) % TRAIL_LEN];
        Uint8 alpha = (Uint8)(30 + 200.0 * (1.0 - (double)k / trail_count));
        SDL_SetRenderDrawColor(rend, 0xFF, 0xFF, 0xFF, alpha);
        SDL_RenderDrawLineF(rend, a.x, a.y, b.x, b.y);
    }
    // Ring around the tracer
    SDL_SetRenderDrawColor(rend, 0xFF, 0xFF, 0xFF, 220);
    double R = p[tracer_i].radius + 4;
    SDL_FPoint ring[25];
    for (int k = 0; k <= 24; k++) {
        double a = 2.0 * M_PI * k / 24;
        ring[k].x = (float)(p[tracer_i].x + R * cos(a));
        ring[k].y = (float)(p[tracer_i].y + R * sin(a));
    }
    SDL_RenderDrawLinesF(rend, ring, 25);
    SDL_SetRenderDrawBlendMode(rend, SDL_BLENDMODE_NONE);
}

static void draw_vessel(SDL_Renderer *rend) {
    SDL_SetRenderDrawColor(rend, 0x55, 0x77, 0xAA, 0xFF);
    for (int xs = 0; xs < cfg.sim_w - 1; xs++) {
        int y0 = (int)para_y(xs), y1 = (int)para_y(xs + 1);
        for (int t = -1; t <= 1; t++)
            SDL_RenderDrawLine(rend, xs, y0 + t, xs + 1, y1 + t);
    }
    {
        int cy2 = (int)ceil_y(sim_time);
        SDL_SetRenderDrawColor(rend, 0x1A, 0x22, 0x3A, 0xFF);
        SDL_Rect piston = {0, 0, cfg.sim_w, cy2};
        SDL_RenderFillRect(rend, &piston);
        SDL_SetRenderDrawColor(rend, 0x55, 0x77, 0xAA, 0xFF);
        for (int t = -1; t <= 1; t++)
            SDL_RenderDrawLine(rend, 0, cy2 + t, cfg.sim_w, cy2 + t);
    }
    double rim_y = (cfg.para_a > 1e-9) ? para_y(cfg.radius) : CFG_PARA_BASE;
    SDL_Rect lw = {0, 0, 2, (int)rim_y}, rw = {cfg.sim_w - 2, 0, 2, (int)rim_y};
    SDL_RenderFillRect(rend, &lw);
    SDL_RenderFillRect(rend, &rw);
}

// ── Tiny 3×5 pixel font ───────────────────────────────────────────────────────
static const Uint8 F[][5] = {
    {0x7,0x5,0x5,0x5,0x7},{0x2,0x6,0x2,0x2,0x7},{0x7,0x1,0x7,0x4,0x7},
    {0x7,0x1,0x7,0x1,0x7},{0x5,0x5,0x7,0x1,0x1},{0x7,0x4,0x7,0x1,0x7},
    {0x7,0x4,0x7,0x5,0x7},{0x7,0x1,0x1,0x1,0x1},{0x7,0x5,0x7,0x5,0x7},
    {0x7,0x5,0x7,0x1,0x7},
    {0x2,0x5,0x7,0x5,0x5},{0x6,0x5,0x6,0x5,0x6},{0x7,0x4,0x4,0x4,0x7},
    {0x6,0x5,0x5,0x5,0x6},{0x7,0x4,0x6,0x4,0x7},{0x7,0x4,0x6,0x4,0x4},
    {0x7,0x4,0x5,0x5,0x7},{0x5,0x5,0x7,0x5,0x5},{0x7,0x2,0x2,0x2,0x7},
    {0x1,0x1,0x1,0x5,0x7},{0x5,0x5,0x6,0x5,0x5},{0x4,0x4,0x4,0x4,0x7},
    {0x5,0x7,0x7,0x5,0x5},{0x5,0x7,0x5,0x5,0x5},{0x2,0x5,0x5,0x5,0x2},
    {0x6,0x5,0x6,0x4,0x4},{0x2,0x5,0x5,0x7,0x3},{0x6,0x5,0x6,0x5,0x5},
    {0x7,0x4,0x7,0x1,0x7},{0x7,0x2,0x2,0x2,0x2},{0x5,0x5,0x5,0x5,0x7},
    {0x5,0x5,0x5,0x5,0x2},{0x5,0x5,0x7,0x7,0x5},{0x5,0x5,0x2,0x5,0x5},
    {0x5,0x5,0x2,0x2,0x2},{0x7,0x1,0x2,0x4,0x7},
    {0x0,0x0,0x0,0x0,0x2},  // 36 '.'
    {0x0,0x0,0x7,0x0,0x0},  // 37 '-'
    {0x1,0x1,0x2,0x4,0x4},  // 38 '/'
    {0x3,0x4,0x4,0x4,0x3},  // 39 '('
    {0x6,0x1,0x1,0x1,0x6},  // 40 ')'
    {0x0,0x0,0x0,0x0,0x0},  // 41 ' '
    {0x0,0x2,0x7,0x2,0x0},  // 42 '+'
    {0x0,0x2,0x0,0x2,0x0},  // 43 ':'
    {0x0,0x7,0x0,0x7,0x0},  // 44 '='
};
static int gl(char c) {
    if (c >= '0' && c <= '9') return c - '0';
    if (c >= 'A' && c <= 'Z') return 10 + (c - 'A');
    if (c >= 'a' && c <= 'z') return 10 + (c - 'a');
    switch (c) {
        case '.': return 36; case '-': return 37; case '/': return 38;
        case '(': return 39; case ')': return 40; case '+': return 42;
        case ':': return 43; case '=': return 44;
        default:  return 41;
    }
}
static void ds(SDL_Renderer *rend, const char *s, int x, int y, int sc,
               Uint8 r, Uint8 g, Uint8 b) {
    SDL_SetRenderDrawColor(rend, r, g, b, 255);
    for (int ci = 0; s[ci]; ci++)
        for (int row = 0; row < 5; row++)
            for (int col = 0; col < 3; col++)
                if (F[gl(s[ci])][row] & (0x4 >> col)) {
                    SDL_Rect px = {x + ci * (3 * sc + sc) + col * sc, y + row * sc, sc, sc};
                    SDL_RenderFillRect(rend, &px);
                }
}
// Right-aligned text at scale 2 (8 px advance)
static void ds_right(SDL_Renderer *rend, const char *s, int xr, int y,
                     Uint8 r, Uint8 g, Uint8 b) {
    ds(rend, s, xr - (int)strlen(s) * 8, y, 2, r, g, b);
}
// Label left, value right, 14 px row
static int row(SDL_Renderer *rend, int x, int xr, int y, const char *label,
               const char *val, Uint8 r, Uint8 g, Uint8 b) {
    ds(rend, label, x, y, 2, 0x88, 0x99, 0xBB);
    ds_right(rend, val, xr, y, r, g, b);
    return y + 14;
}

// ── Plot helpers ──────────────────────────────────────────────────────────────
static void plot_box(SDL_Renderer *rend, SDL_Rect b) {
    SDL_SetRenderDrawColor(rend, 0x06, 0x08, 0x10, 0xFF);
    SDL_RenderFillRect(rend, &b);
    SDL_SetRenderDrawBlendMode(rend, SDL_BLENDMODE_BLEND);
    SDL_SetRenderDrawColor(rend, 0x22, 0x33, 0x55, 80);
    for (int f = 1; f <= 3; f++)
        SDL_RenderDrawLine(rend, b.x, b.y + b.h - b.h * f / 4, b.x + b.w, b.y + b.h - b.h * f / 4);
    SDL_SetRenderDrawBlendMode(rend, SDL_BLENDMODE_NONE);
}
static void plot_border(SDL_Renderer *rend, SDL_Rect b) {
    SDL_SetRenderDrawColor(rend, 0x2A, 0x3A, 0x55, 0xFF);
    SDL_RenderDrawRect(rend, &b);
}
static inline int plot_y(SDL_Rect b, double v, double lo, double hi) {
    double t = (hi > lo) ? (v - lo) / (hi - lo) : 0.5;
    int y = b.y + b.h - 1 - (int)(t * (b.h - 2));
    if (y < b.y) y = b.y;
    if (y >= b.y + b.h) y = b.y + b.h - 1;
    return y;
}
static void plot_series(SDL_Renderer *rend, SDL_Rect b, const Series *s,
                        double lo, double hi, Uint8 r, Uint8 g, Uint8 bl) {
    if (s->count < 2) return;
    int cols = s->count < b.w ? s->count : b.w;
    SDL_SetRenderDrawColor(rend, r, g, bl, 255);
    int px = -1, py = -1;
    for (int c = 0; c < cols; c++) {
        int x = b.x + b.w - 1 - c, y = plot_y(b, series_at(s, c), lo, hi);
        if (px >= 0) SDL_RenderDrawLine(rend, px, py, x, y);
        px = x; py = y;
    }
}
static void plot_hline(SDL_Renderer *rend, SDL_Rect b, double v, double lo, double hi,
                       Uint8 r, Uint8 g, Uint8 bl, Uint8 a) {
    int y = plot_y(b, v, lo, hi);
    SDL_SetRenderDrawBlendMode(rend, SDL_BLENDMODE_BLEND);
    SDL_SetRenderDrawColor(rend, r, g, bl, a);
    SDL_RenderDrawLine(rend, b.x, y, b.x + b.w - 1, y);
    SDL_SetRenderDrawBlendMode(rend, SDL_BLENDMODE_NONE);
}
static void pad_range(double *lo, double *hi, double frac) {
    double span = *hi - *lo;
    if (span < 1e-9) span = fabs(*hi) > 1e-9 ? fabs(*hi) * 0.02 : 1.0;
    *lo -= span * frac; *hi += span * frac;
}

// ── Slider system ─────────────────────────────────────────────────────────────
static int g_slider_panel_top = 500;

typedef struct {
    const char *label;
    double min, max;
    double *dvalue;
    int    *ivalue;
    int     needs_reset;
} Slider;

static Slider sliders[10];
static int    n_sliders = 0;
static int    active_slider = -1;

static void sliders_init(void) {
    n_sliders = 0;
#define SD(lbl,mn,mx,fld,rst) sliders[n_sliders++]=(Slider){lbl,mn,mx,&cfg.fld,NULL,rst}
#define SI(lbl,mn,mx,fld,rst) sliders[n_sliders++]=(Slider){lbl,mn,mx,NULL,&cfg.fld,rst}
    SI("PARTICLES",   1, 5000, n_particles, 1);
    SD("GRAVITY",     0, 2000, gravity,     0);
    SD("RADIUS",      1,   20, radius,      1);
    SD("RESTITUTION", 0,    1, restitution, 0);
    SD("MONO SPEED", 10, 1000, mono_speed,  1);
    SD("PARA A",      0, 0.01, para_a,      1);
    SD("PISTON AMP",  0,  150, ceil_amp,    0);
    SD("PISTON FREQ", 0,   10, ceil_freq,   0);
    SD("HEAVY FRAC",  0,    1, heavy_frac,  1);
    SD("MASS RATIO",  1,   20, mass_ratio,  1);
#undef SD
#undef SI
}
static double slider_get(int i) {
    return sliders[i].dvalue ? *sliders[i].dvalue : (double)*sliders[i].ivalue;
}
static void slider_put(int i, double v) {
    if (sliders[i].dvalue) *sliders[i].dvalue = v;
    else *sliders[i].ivalue = (int)round(v);
}
static SDL_Rect slider_track(int i, int top) {
    return (SDL_Rect){cfg.sim_w + 8, top + 18 + i * 30 + 16, PANEL_W - 16, 5};
}
static void draw_sliders(SDL_Renderer *rend, int top) {
    ds(rend, "PARAMETERS", cfg.sim_w + 8, top, 2, 0xAA, 0xBB, 0xCC);
    for (int i = 0; i < n_sliders; i++) {
        Slider *s = &sliders[i];
        SDL_Rect tr = slider_track(i, top);
        int ly = tr.y - 13;
        ds(rend, s->label, tr.x, ly, 2, 0x88, 0x99, 0xBB);
        char buf[16];
        double v = slider_get(i);
        if (s->ivalue)         snprintf(buf, sizeof(buf), "%d", (int)v);
        else if (s->max <= 0.1)  snprintf(buf, sizeof(buf), "%.4f", v);
        else if (s->max <= 1.0)  snprintf(buf, sizeof(buf), "%.3f", v);
        else if (s->max <= 10.0) snprintf(buf, sizeof(buf), "%.2f", v);
        else                     snprintf(buf, sizeof(buf), "%.1f", v);
        ds_right(rend, buf, tr.x + tr.w, ly, 0xFF, 0xDD, 0x88);
        SDL_SetRenderDrawColor(rend, 0x22, 0x2A, 0x3A, 0xFF);
        SDL_RenderFillRect(rend, &tr);
        double frac = (v - s->min) / (s->max - s->min);
        if (frac < 0) frac = 0;
        if (frac > 1) frac = 1;
        SDL_Rect fill = tr; fill.w = (int)(frac * tr.w);
        SDL_SetRenderDrawColor(rend,
            i == active_slider ? 0x88 : 0x44,
            i == active_slider ? 0xCC : 0x88,
            i == active_slider ? 0xFF : 0xCC, 0xFF);
        if (fill.w > 0) SDL_RenderFillRect(rend, &fill);
        SDL_Rect thumb = {tr.x + fill.w - 4, tr.y - 4, 8, 13};
        SDL_SetRenderDrawColor(rend, 0xDD, 0xEE, 0xFF, 0xFF);
        SDL_RenderFillRect(rend, &thumb);
    }
}
static int slider_hit(int mx, int my, int top) {
    for (int i = 0; i < n_sliders; i++) {
        SDL_Rect tr = slider_track(i, top);
        SDL_Rect h  = {tr.x, tr.y - 6, tr.w, tr.h + 12};
        if (mx >= h.x && mx < h.x + h.w && my >= h.y && my < h.y + h.h) return i;
    }
    return -1;
}
static int slider_set(int idx, int mx, int top) {
    Slider *s = &sliders[idx];
    SDL_Rect tr = slider_track(idx, top);
    double frac = (double)(mx - tr.x) / tr.w;
    if (frac < 0) frac = 0;
    if (frac > 1) frac = 1;
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

// ── Stats panel (middle column) ───────────────────────────────────────────────
static void draw_panel(SDL_Renderer *rend, int n, double total_e, double ke, double pe) {
    SDL_SetRenderDrawColor(rend, 0x0C, 0x0F, 0x1C, 0xFF);
    SDL_Rect bg = {cfg.sim_w, 0, PANEL_W, cfg.sim_h};
    SDL_RenderFillRect(rend, &bg);
    SDL_SetRenderDrawColor(rend, 0x44, 0x55, 0x77, 0xFF);
    SDL_RenderDrawLine(rend, cfg.sim_w, 0, cfg.sim_w, cfg.sim_h);

    int px = cfg.sim_w + 8, ty = 10;
    char buf[48];

    // REC indicator
    {
        int rx = cfg.sim_w + PANEL_W - 8;
        if (rec_active) {
            if (((SDL_GetTicks() - rec_blink) % 900) < 600) {
                SDL_SetRenderDrawColor(rend, 0xEE, 0x22, 0x22, 0xFF);
                for (int dy = -5; dy <= 5; dy++) for (int dx = -5; dx <= 5; dx++)
                    if (dx * dx + dy * dy <= 25) SDL_RenderDrawPoint(rend, rx - 6 + dx, ty + 5 + dy);
            }
            ds(rend, "REC", rx - 38, ty, 2, 0xEE, 0x33, 0x33);
            snprintf(buf, sizeof(buf), "%ldf", rec_frames);
            ds_right(rend, buf, rx, ty + 14, 0xAA, 0x44, 0x44);
        } else {
            SDL_SetRenderDrawColor(rend, 0x33, 0x33, 0x33, 0xFF);
            for (int dy = -4; dy <= 4; dy++) for (int dx = -4; dx <= 4; dx++)
                if (dx * dx + dy * dy <= 16) SDL_RenderDrawPoint(rend, rx - 6 + dx, ty + 4 + dy);
            ds(rend, "V-REC", rx - 46, ty, 2, 0x44, 0x44, 0x55);
        }
    }

    ds(rend, "V2.1", px, ty, 2, 0x33, 0x88, 0xFF);
    ty += 18;

    snprintf(buf, sizeof(buf), "%dx%d", cfg.sim_w, cfg.sim_h);
    ty = draw_stat(rend, px, ty, "BOX SIZE (cfg sim-w/h)", buf, 0x88, 0x99, 0xAA) + 4;
    snprintf(buf, sizeof(buf), "%d", n);
    ty = draw_stat(rend, px, ty, "PARTICLES", buf, 0xFF, 0xFF, 0xFF) + 4;
    ty = row(rend, px, cfg.sim_w + PANEL_W - 8, ty, "COLOR MODE (C)", color_mode_names[color_mode], 0xFF, 0xDD, 0xFF) + 4;
    snprintf(buf, sizeof(buf), "%.4f", (eng_e0 > 0) ? total_e / eng_e0 : 1.0);
    ty = draw_stat(rend, px, ty, "TOTAL ENERGY (E/E0)", buf, 0xFF, 0xCC, 0x44) + 4;
    snprintf(buf, sizeof(buf), "%.4f", (total_e != 0) ? ke / total_e : 0.0);
    ty = draw_stat(rend, px, ty, "KINETIC (KE/E)", buf, 0x44, 0xFF, 0xAA) + 4;
    snprintf(buf, sizeof(buf), "%.4f", (total_e != 0) ? pe / total_e : 0.0);
    ty = draw_stat(rend, px, ty, "POTENTIAL (PE/E)", buf, 0x88, 0xAA, 0xFF) + 4;
    snprintf(buf, sizeof(buf), "%.1f", (n > 0) ? ke / n : 0.0);
    ty = draw_stat(rend, px, ty, "TEMPERATURE (KE/N)", buf, 0xFF, 0x66, 0x44) + 4;
    snprintf(buf, sizeof(buf), "%.1f", species_n[0] > 0 ? species_ke[0] / species_n[0] : 0.0);
    ty = draw_stat(rend, px, ty, "LIGHT KE/N", buf, 0x33, 0xCC, 0xFF) + 4;
    snprintf(buf, sizeof(buf), "%.1f", species_n[1] > 0 ? species_ke[1] / species_n[1] : 0.0);
    ty = draw_stat(rend, px, ty, "HEAVY KE/N", buf, 0xFF, 0x88, 0x33) + 4;
    snprintf(buf, sizeof(buf), "%ld", collision_count);
    ty = draw_stat(rend, px, ty, "COLLISIONS", buf, 0xCC, 0x88, 0xFF) + 10;

    // Energy trace
    ds(rend, "ENERGY / E0", px, ty, 2, 0xAA, 0xBB, 0xCC);
    ty += 13;
    SDL_Rect gb = {cfg.sim_w + 4, ty, PANEL_W - 8, 70};
    plot_box(rend, gb);
    double emin, emax;
    series_minmax(&s_energy, &emin, &emax);
    if (emax < 1.0) emax = 1.0;
    if (emin > 1.0) emin = 1.0;
    if (emax - emin < 0.002) { emax = 1.001; emin = 0.999; }
    plot_hline(rend, gb, 1.0, emin, emax, 0xFF, 0xFF, 0xFF, 60);
    if (s_energy.count > 1) {
        int cols = s_energy.count < gb.w ? s_energy.count : gb.w;
        double range = emax - emin;
        for (int c = 0; c < cols; c++) {
            double v = series_at(&s_energy, c);
            int by = plot_y(gb, v, emin, emax), bx = gb.x + gb.w - 1 - c;
            double dev = fabs(v - 1.0) / (range * 0.5 + 1e-9);
            if (dev > 1) dev = 1;
            SDL_SetRenderDrawColor(rend, (Uint8)(dev * 255), (Uint8)((1 - dev) * 200), 0x18, 0xFF);
            SDL_RenderDrawLine(rend, bx, by, bx, gb.y + gb.h - 1);
        }
    }
    plot_border(rend, gb);
    if (s_energy.count > 0) {
        double pct = (series_at(&s_energy, 0) - 1.0) * 100.0;
        snprintf(buf, sizeof(buf), "%+.3f PCT", pct);
        ds(rend, "DRIFT", px, ty + gb.h + 2, 2, 0xAA, 0xBB, 0xCC);
        Uint8 dcr = (Uint8)(fabs(pct) > 1.0 ? 220 : fabs(pct) > 0.1 ? 180 : 100);
        Uint8 dcg = (Uint8)(fabs(pct) > 1.0 ? 80 : fabs(pct) > 0.1 ? 180 : 220);
        ds(rend, buf, px + 48, ty + gb.h + 2, 2, dcr, dcg, 0x44);
    }
    ty += gb.h + 18;

    // KE histogram
    ds(rend, "KE DISTRIBUTION vs MB", px, ty, 2, 0xAA, 0xBB, 0xCC);
    ty += 13;
    SDL_Rect hb = {cfg.sim_w + 4, ty, PANEL_W - 8, 85};
    plot_box(rend, hb);
    double bin_w = (double)hb.w / HIST_BINS;
    for (int b = 0; b < HIST_BINS; b++) {
        int bx = (int)(hb.x + b * bin_w), bw = (int)bin_w - 1;
        if (bw < 1) bw = 1;
        int bar_h = (int)(hist_smooth[b] * (hb.h - 2));
        if (bar_h > 0) {
            SDL_SetRenderDrawBlendMode(rend, SDL_BLENDMODE_BLEND);
            SDL_SetRenderDrawColor(rend, 0x22, 0xAA, 0xCC, 180);
            SDL_Rect br = {bx, hb.y + hb.h - bar_h, bw, bar_h};
            SDL_RenderFillRect(rend, &br);
            SDL_SetRenderDrawBlendMode(rend, SDL_BLENDMODE_NONE);
        }
        if (b > 0) {
            int y0 = (int)(hb.y + hb.h - hist_mb[b - 1] * (hb.h - 2));
            int y1 = (int)(hb.y + hb.h - hist_mb[b] * (hb.h - 2));
            SDL_SetRenderDrawColor(rend, 0xFF, 0xDD, 0x33, 0xFF);
            SDL_RenderDrawLine(rend, (int)(hb.x + (b - 0.5) * bin_w), y0,
                                     (int)(hb.x + (b + 0.5) * bin_w), y1);
        }
    }
    plot_border(rend, hb);
    ty += hb.h + 6;

    SDL_SetRenderDrawColor(rend, 0x22, 0xAA, 0xCC, 0xFF);
    SDL_Rect lc = {px, ty + 2, 8, 8}; SDL_RenderFillRect(rend, &lc);
    ds(rend, "MEASURED", px + 11, ty, 2, 0x44, 0xCC, 0xEE);
    ty += 14;
    SDL_SetRenderDrawColor(rend, 0xFF, 0xDD, 0x33, 0xFF);
    SDL_RenderDrawLine(rend, px, ty + 4, px + 8, ty + 4);
    ds(rend, "MAXWELL-BOLTZMANN", px + 11, ty, 2, 0xDD, 0xBB, 0x44);
    ty += 20;
    SDL_SetRenderDrawColor(rend, 0x33, 0xCC, 0xFF, 0xFF);
    SDL_Rect sl = {px, ty + 2, 8, 8}; SDL_RenderFillRect(rend, &sl);
    ds(rend, "LIGHT", px + 11, ty, 2, 0x33, 0xCC, 0xFF);
    SDL_SetRenderDrawColor(rend, 0xFF, 0x88, 0x33, 0xFF);
    SDL_Rect sh = {px + 60, ty + 2, 8, 8}; SDL_RenderFillRect(rend, &sh);
    ds(rend, "HEAVY", px + 71, ty, 2, 0xFF, 0x88, 0x33);
    ty += 14;

    SDL_SetRenderDrawColor(rend, 0x2A, 0x3A, 0x55, 0xFF);
    SDL_RenderDrawLine(rend, cfg.sim_w + 4, ty, cfg.sim_w + PANEL_W - 4, ty);
    ty += 6;

    g_slider_panel_top = ty;
    draw_sliders(rend, ty);

    ds(rend, "R-RAND M-MONO C-COLOR T-TRACER V-REC Q-QUIT",
       px, cfg.sim_h - 14, 1, 0x55, 0x66, 0x88);
}

// ── Analysis panel (right column) ─────────────────────────────────────────────
static void draw_analysis(SDL_Renderer *rend, int n, double ke) {
    int ax = cfg.sim_w + PANEL_W;
    int px = ax + 8, xr = ax + ANALYSIS_W - 8, aw = ANALYSIS_W - 16;
    int ty = 10;
    char buf[48];
    double kT = (n > 0) ? ke / n : 0.0;

    SDL_SetRenderDrawColor(rend, 0x0A, 0x0D, 0x18, 0xFF);
    SDL_Rect bg = {ax, 0, ANALYSIS_W, cfg.sim_h};
    SDL_RenderFillRect(rend, &bg);
    SDL_SetRenderDrawColor(rend, 0x44, 0x55, 0x77, 0xFF);
    SDL_RenderDrawLine(rend, ax, 0, ax, cfg.sim_h);

    ds(rend, "ANALYSIS", px, ty, 2, 0x33, 0x88, 0xFF);
    ty += 16;

    // ── 1. Equipartition ────────────────────────────────────────────────────
    ds(rend, "EQUIPARTITION  KE/N", px, ty, 2, 0xAA, 0xBB, 0xCC);
    ty += 13;
    {
        SDL_Rect b = {px, ty, aw, 80};
        plot_box(rend, b);
        double lo0, hi0, lo1, hi1;
        series_minmax(&s_eq[0], &lo0, &hi0);
        series_minmax(&s_eq[1], &lo1, &hi1);
        double hi = fmax(fmax(hi0, hi1), kT) * 1.1, lo = 0.0;
        if (hi < 1) hi = 1;
        plot_hline(rend, b, kT, lo, hi, 0xFF, 0xFF, 0xFF, 70);
        plot_series(rend, b, &s_eq[0], lo, hi, 0x33, 0xCC, 0xFF);
        plot_series(rend, b, &s_eq[1], lo, hi, 0xFF, 0x88, 0x33);
        plot_border(rend, b);
        ty += b.h + 4;
        snprintf(buf, sizeof(buf), "L %.0f  H %.0f  MEAN %.0f",
                 species_n[0] > 0 ? species_ke[0] / species_n[0] : 0.0,
                 species_n[1] > 0 ? species_ke[1] / species_n[1] : 0.0, kT);
        ds(rend, buf, px, ty, 2, 0x99, 0xAA, 0xBB);
        ty += 16;
    }

    // ── 2. H-function ───────────────────────────────────────────────────────
    ds(rend, "H-FUNCTION (KE MARGINAL)", px, ty, 2, 0xAA, 0xBB, 0xCC);
    ty += 13;
    {
        SDL_Rect b = {px, ty, aw, 60};
        plot_box(rend, b);
        double lo, hi;
        series_minmax(&s_h, &lo, &hi);
        if (h_eq < lo) lo = h_eq;
        if (h_eq > hi) hi = h_eq;
        pad_range(&lo, &hi, 0.1);
        plot_hline(rend, b, h_eq, lo, hi, 0xFF, 0xDD, 0x33, 120);
        plot_series(rend, b, &s_h, lo, hi, 0xCC, 0x88, 0xFF);
        plot_border(rend, b);
        ty += b.h + 4;
        snprintf(buf, sizeof(buf), "H %.3f   H-EQ %.3f",
                 s_h.count ? series_at(&s_h, 0) : 0.0, h_eq);
        ds(rend, buf, px, ty, 2, 0x99, 0xAA, 0xBB);
        ty += 16;
    }

    // ── 3. Barometric density profile ──────────────────────────────────────
    ds(rend, "DENSITY N(Y) vs EXP(-MGH/KT)", px, ty, 2, 0xAA, 0xBB, 0xCC);
    ty += 13;
    {
        SDL_Rect b = {px, ty, aw, 150};
        SDL_SetRenderDrawColor(rend, 0x06, 0x08, 0x10, 0xFF);
        SDL_RenderFillRect(rend, &b);
        int cxm = b.x + b.w / 2, half = b.w / 2 - 4;
        double rowh = (double)b.h / NB_Y;
        double y_top = ceil_y(sim_time), y_bot = CFG_PARA_BASE;
        double span = y_bot - y_top; if (span < 1) span = 1;
        for (int s = 0; s < 2; s++) {
            // Theory: n_b ∝ exp(-m g h_b / kT) · width(y_b), normalised to N_s
            double th[NB_Y], tsum = 0, msum = 0, mass = (s == 0) ? 1.0 : cfg.mass_ratio;
            for (int k = 0; k < NB_Y; k++) {
                double yb = y_top + (k + 0.5) / NB_Y * span, h = y_bot - yb;
                double w = vessel_width(yb);
                th[k] = (kT > 1e-9) ? exp(-mass * cfg.gravity * h / kT) * w : w;
                tsum += th[k];
                msum += dens_meas[s][k];
            }
            double mx = 0;
            for (int k = 0; k < NB_Y; k++) {
                th[k] = (tsum > 0) ? th[k] / tsum * species_n[s] : 0;
                if (th[k] > mx) mx = th[k];
                if (dens_meas[s][k] > mx) mx = dens_meas[s][k];
            }
            if (mx < 1e-9) mx = 1;
            int dir = (s == 0) ? -1 : 1;
            Uint8 cr = s ? 0xFF : 0x33, cg = s ? 0x88 : 0xCC, cb = s ? 0x33 : 0xFF;
            SDL_SetRenderDrawBlendMode(rend, SDL_BLENDMODE_BLEND);
            SDL_SetRenderDrawColor(rend, cr, cg, cb, 150);
            for (int k = 0; k < NB_Y; k++) {
                int len = (int)(dens_meas[s][k] / mx * half);
                int y0 = b.y + (int)(k * rowh), hh = (int)rowh - 1;
                if (hh < 1) hh = 1;
                SDL_Rect r = {dir < 0 ? cxm - len : cxm + 1, y0, len, hh};
                if (len > 0) SDL_RenderFillRect(rend, &r);
            }
            SDL_SetRenderDrawBlendMode(rend, SDL_BLENDMODE_NONE);
            SDL_SetRenderDrawColor(rend, 0xFF, 0xDD, 0x33, 0xFF);
            int pxl = -1, pyl = -1;
            for (int k = 0; k < NB_Y; k++) {
                int x = cxm + dir * (int)(th[k] / mx * half);
                int y = b.y + (int)((k + 0.5) * rowh);
                if (pxl >= 0) SDL_RenderDrawLine(rend, pxl, pyl, x, y);
                pxl = x; pyl = y;
            }
            (void)msum;
        }
        SDL_SetRenderDrawColor(rend, 0x2A, 0x3A, 0x55, 0xFF);
        SDL_RenderDrawLine(rend, cxm, b.y, cxm, b.y + b.h);
        plot_border(rend, b);
        ty += b.h + 4;
        double hs_l = (cfg.gravity > 0 && kT > 0) ? kT / cfg.gravity : 0;
        double hs_h = (cfg.gravity > 0 && kT > 0) ? kT / (cfg.mass_ratio * cfg.gravity) : 0;
        ds(rend, "LIGHT", px, ty, 2, 0x33, 0xCC, 0xFF);
        ds(rend, "(TOP = UP)  MB THEORY", px + 60, ty, 2, 0x99, 0xAA, 0xBB);
        ds_right(rend, "HEAVY", xr, ty, 0xFF, 0x88, 0x33);
        ty += 14;
        snprintf(buf, sizeof(buf), "SCALE HEIGHT KT/MG  L %.0f  H %.0f", hs_l, hs_h);
        ds(rend, buf, px, ty, 2, 0x99, 0xAA, 0xBB);
        ty += 16;
    }

    // ── 4. Pressure ─────────────────────────────────────────────────────────
    ds(rend, "WALL PRESSURE (2D: P=F/L)", px, ty, 2, 0xAA, 0xBB, 0xCC);
    ty += 13;
    snprintf(buf, sizeof(buf), "%.2f", P_smooth);
    ty = row(rend, px, xr, ty, "P MEASURED", buf, 0xFF, 0xDD, 0x88);
    snprintf(buf, sizeof(buf), "%.2f", P_ideal);
    ty = row(rend, px, xr, ty, "P IDEAL  NKT/A", buf, 0xFF, 0xDD, 0x88);
    // Z only means something without gravity (the floor also carries the weight)
    bool zok = cfg.gravity < 1e-9;
    snprintf(buf, sizeof(buf), "%.3f", Z_meas);
    ty = row(rend, px, xr, ty, zok ? "Z = PA/NKT" : "Z = PA/NKT (NEEDS G=0)",
             buf, zok ? 0x44 : 0x33, zok ? 0xFF : 0x44, zok ? 0xAA : 0x44);
    snprintf(buf, sizeof(buf), "%.3f  ETA %.3f", Z_hend, packing);
    ty = row(rend, px, xr, ty, "Z HENDERSON", buf, zok ? 0xDD : 0x55, zok ? 0xBB : 0x55, 0x44);
    ty += 4;

    // ── 5. P-V diagram ──────────────────────────────────────────────────────
    ds(rend, "P-V LOOP (PISTON)", px, ty, 2, 0xAA, 0xBB, 0xCC);
    ty += 13;
    {
        SDL_Rect b = {px, ty, aw, 90};
        plot_box(rend, b);
        if (s_P.count > 1) {
            double plo, phi, vlo, vhi;
            series_minmax(&s_P, &plo, &phi);
            series_minmax(&s_V, &vlo, &vhi);
            pad_range(&plo, &phi, 0.1);
            pad_range(&vlo, &vhi, 0.1);
            int cols = s_P.count < 160 ? s_P.count : 160;
            SDL_SetRenderDrawBlendMode(rend, SDL_BLENDMODE_BLEND);
            int pxl = -1, pyl = -1;
            for (int c = cols - 1; c >= 0; c--) {
                double P = series_at(&s_P, c), V = series_at(&s_V, c);
                int x = b.x + 1 + (int)((V - vlo) / (vhi - vlo) * (b.w - 3));
                int y = plot_y(b, P, plo, phi);
                Uint8 a = (Uint8)(40 + 200.0 * (1.0 - (double)c / cols));
                SDL_SetRenderDrawColor(rend, 0x44, 0xFF, 0xAA, a);
                if (pxl >= 0) SDL_RenderDrawLine(rend, pxl, pyl, x, y);
                pxl = x; pyl = y;
            }
            SDL_SetRenderDrawBlendMode(rend, SDL_BLENDMODE_NONE);
            SDL_SetRenderDrawColor(rend, 0xFF, 0xFF, 0xFF, 0xFF);
            SDL_Rect dot = {pxl - 2, pyl - 2, 5, 5};
            SDL_RenderFillRect(rend, &dot);
        }
        plot_border(rend, b);
        ty += b.h + 4;
        if (cfg.ceil_amp > 1e-9) {
            double pct = (eng_e0 > 0) ? 100.0 * piston_work_last / eng_e0 : 0.0;
            if (fabs(piston_work_last) >= 1e5)
                snprintf(buf, sizeof(buf), "W/CYC %+.2e (%+.2f PCT E0) %ld", piston_work_last, pct, piston_cycles);
            else
                snprintf(buf, sizeof(buf), "W/CYC %+.1f (%+.2f PCT E0) %ld", piston_work_last, pct, piston_cycles);
        } else {
            snprintf(buf, sizeof(buf), "PISTON OFF - SET PISTON AMP");
        }
        ds(rend, buf, px, ty, 2, 0x99, 0xAA, 0xBB);
        ty += 16;
    }

    // ── 6. Tracer MSD ───────────────────────────────────────────────────────
    ds(rend, "HEAVY-PARTICLE MSD  (T RESET)", px, ty, 2, 0xAA, 0xBB, 0xCC);
    ty += 13;
    {
        SDL_Rect b = {px, ty, aw, 70};
        plot_box(rend, b);
        double lo = 0, hi, tmp;
        series_minmax(&s_msd, &tmp, &hi);
        hi *= 1.1; if (hi < 1) hi = 1;
        plot_series(rend, b, &s_msd, lo, hi, 0xFF, 0xFF, 0xFF);
        plot_border(rend, b);
        ty += b.h + 4;
        double msd = s_msd.count ? series_at(&s_msd, 0) : 0.0;
        snprintf(buf, sizeof(buf), "MSD %.0f PX2   RMS %.1f PX", msd, sqrt(msd));
        ds(rend, buf, px, ty, 2, 0x99, 0xAA, 0xBB);
        ty += 16;
    }

    // ── 7. Mean free path ───────────────────────────────────────────────────
    ds(rend, "KINETIC THEORY", px, ty, 2, 0xAA, 0xBB, 0xCC);
    ty += 13;
    snprintf(buf, sizeof(buf), "%.1f /S", coll_rate);
    ty = row(rend, px, xr, ty, "COLL RATE / PARTICLE", buf, 0xCC, 0x88, 0xFF);
    snprintf(buf, sizeof(buf), "%.1f PX/S", mean_speed);
    ty = row(rend, px, xr, ty, "MEAN SPEED", buf, 0xCC, 0x88, 0xFF);
    snprintf(buf, sizeof(buf), "%.1f PX", mfp_meas);
    ty = row(rend, px, xr, ty, "MFP MEASURED", buf, 0x44, 0xFF, 0xAA);
    snprintf(buf, sizeof(buf), "%.1f PX", mfp_theory);
    ty = row(rend, px, xr, ty, "MFP ENSKOG (LOCAL N)", buf, 0xDD, 0xBB, 0x44);
    ty += 4;

    // ── 8. Performance ──────────────────────────────────────────────────────
    ds(rend, "PERFORMANCE", px, ty, 2, 0xAA, 0xBB, 0xCC);
    ty += 13;
    snprintf(buf, sizeof(buf), "%.2f MS", perf_step_ms);
    ty = row(rend, px, xr, ty, "PHYSICS / FRAME", buf, 0x99, 0xAA, 0xBB);
    int threads = 1;
#ifdef _OPENMP
    threads = omp_get_max_threads();
#endif
    snprintf(buf, sizeof(buf), "%.0f FPS  %d THR  %dx%d", perf_fps, threads, g_substeps, g_passes);
    ty = row(rend, px, xr, ty, "FRAME", buf, 0x99, 0xAA, 0xBB);
}

// ── One full frame ────────────────────────────────────────────────────────────
static void render_frame(SDL_Renderer *rend, int n) {
    SDL_SetRenderDrawColor(rend, 0x0D, 0x0F, 0x1A, 0xFF);
    SDL_RenderClear(rend);

    double max_val = 1.0;
    if (color_mode == 1) {
        for (int i = 0; i < n; i++) {
            double s = sqrt(p[i].vx * p[i].vx + p[i].vy * p[i].vy);
            if (s > max_val) max_val = s;
        }
    } else if (color_mode == 2) {
        for (int i = 0; i < n; i++) {
            double k = 0.5 * p[i].mass * (p[i].vx * p[i].vx + p[i].vy * p[i].vy);
            if (k > max_val) max_val = k;
        }
    }
    draw_particles(rend, n, max_val);
    draw_tracer(rend, n);
    draw_vessel(rend);

    double ke, pe;
    double te = compute_energy(n, &ke, &pe);
    draw_panel(rend, n, te, ke, pe);
    draw_analysis(rend, n, ke);
}

static void trail_push(int n) {
    if (tracer_i >= n) return;
    trail[trail_head] = (SDL_FPoint){(float)p[tracer_i].x, (float)p[tracer_i].y};
    trail_head = (trail_head + 1) % TRAIL_LEN;
    if (trail_count < TRAIL_LEN) trail_count++;
}

// ── Config file search ────────────────────────────────────────────────────────
// Without --config, the first readable file in this order is used:
//   1. ./ballsim_v2.cfg                              (project-local override)
//   2. $XDG_CONFIG_HOME/ballsim/ballsim_v2.cfg       (~/.config — user edits)
//   3. $XDG_DATA_HOME/ballsim/ballsim_v2.cfg         (~/.local/share — installed default)
//   4. DATADIR/ballsim_v2.cfg                        (compile-time, from `make install`)
// Matching the XDG convention: `make install` puts the pristine default in
// share/, and the user copies it to ~/.config/ballsim/ to customise.
#ifndef BALLSIM_VERSION
#define BALLSIM_VERSION "2.1"
#endif
#ifndef DATADIR
#define DATADIR ""
#endif
#define CFG_NAME "ballsim_v2.cfg"

static bool file_readable(const char *path) {
    FILE *f = fopen(path, "r");
    if (!f) return false;
    fclose(f);
    return true;
}

static const char *find_config(char *buf, size_t sz) {
    if (file_readable(CFG_NAME)) return CFG_NAME;
    const char *home = getenv("HOME");
    const char *xc = getenv("XDG_CONFIG_HOME"), *xd = getenv("XDG_DATA_HOME");
    if (xc && *xc) snprintf(buf, sz, "%s/ballsim/" CFG_NAME, xc);
    else if (home) snprintf(buf, sz, "%s/.config/ballsim/" CFG_NAME, home);
    else buf[0] = 0;
    if (buf[0] && file_readable(buf)) return buf;
    if (xd && *xd) snprintf(buf, sz, "%s/ballsim/" CFG_NAME, xd);
    else if (home) snprintf(buf, sz, "%s/.local/share/ballsim/" CFG_NAME, home);
    else buf[0] = 0;
    if (buf[0] && file_readable(buf)) return buf;
    if (DATADIR[0]) {
        snprintf(buf, sz, "%s/" CFG_NAME, DATADIR);
        if (file_readable(buf)) return buf;
    }
    return CFG_NAME;   // nothing found; load_config will report and use defaults
}

// ── Main ──────────────────────────────────────────────────────────────────────
static void usage(void) {
    printf("Usage: ballsim_v2 [options]\n"
           "  --config FILE          load config from FILE\n"
           "  --write-config [FILE]  write default config and exit (default: " CFG_NAME ")\n"
           "  --bench N STEPS        headless benchmark, no window\n"
           "  --screenshot FILE K    run K steps offscreen, save BMP, exit\n"
           "  --substeps K           Verlet sub-steps per 1/120 s (default %d)\n"
           "  --passes P             resolution passes per sub-step (default %d)\n"
           "  --threads T            OpenMP thread count\n"
           "  --seed S               RNG seed (default: time)\n"
           "  --version              print version and exit\n"
           "Config search order without --config: ./" CFG_NAME ", "
           "$XDG_CONFIG_HOME/ballsim/, $XDG_DATA_HOME/ballsim/%s%s\n",
           g_substeps, g_passes, DATADIR[0] ? ", " : "", DATADIR);
}

int main(int argc, char **argv) {
    const char *config_path = NULL;
    char        cfg_buf[1024];
    const char *shot_path = NULL;
    int bench_n = 0, bench_steps = 0, shot_steps = 0, max_frames = 0;
    unsigned seed = (unsigned)time(NULL);
    bool write_cfg = false; const char *write_path = "ballsim_v2.cfg";

    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "--config") && i + 1 < argc) config_path = argv[++i];
        else if (!strcmp(argv[i], "--write-config")) {
            write_cfg = true;
            if (i + 1 < argc && argv[i + 1][0] != '-') write_path = argv[++i];
        }
        else if (!strcmp(argv[i], "--bench") && i + 2 < argc) {
            bench_n = atoi(argv[++i]); bench_steps = atoi(argv[++i]);
        }
        else if (!strcmp(argv[i], "--screenshot") && i + 2 < argc) {
            shot_path = argv[++i]; shot_steps = atoi(argv[++i]);
        }
        else if (!strcmp(argv[i], "--substeps") && i + 1 < argc) g_substeps = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--passes")   && i + 1 < argc) g_passes   = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--seed")     && i + 1 < argc) seed = (unsigned)atoi(argv[++i]);
        else if (!strcmp(argv[i], "--frames")   && i + 1 < argc) max_frames = atoi(argv[++i]); // test: quit after K frames
        else if (!strcmp(argv[i], "--threads")  && i + 1 < argc) {
#ifdef _OPENMP
            omp_set_num_threads(atoi(argv[++i]));
#else
            ++i; fprintf(stderr, "Built without OpenMP; --threads ignored.\n");
#endif
        }
        else if (!strcmp(argv[i], "--help") || !strcmp(argv[i], "-h")) { usage(); return 0; }
        else if (!strcmp(argv[i], "--version")) { printf("ballsim_v2 " BALLSIM_VERSION "\n"); return 0; }
        else { fprintf(stderr, "Unknown option: %s\n", argv[i]); usage(); return 1; }
    }
    if (g_substeps < 1) g_substeps = 1;
    if (g_passes < 1)   g_passes   = 1;
    srand(seed);

    if (write_cfg) { config_derive(); write_default_config(write_path); return 0; }
    if (!config_path) config_path = find_config(cfg_buf, sizeof(cfg_buf));
    load_config(config_path);
    config_derive();

    if (bench_n > 0 && bench_steps > 0) return run_bench(bench_n, bench_steps);

    cell_list_init();
#ifdef _OPENMP
    fprintf(stderr, "OpenMP: %d threads.\n", omp_get_max_threads());
#endif

    if (shot_path || max_frames > 0) SDL_SetHint("SDL_VIDEODRIVER", "dummy");
    SDL_SetHint(SDL_HINT_RENDER_SCALE_QUALITY, "1");
    if (SDL_Init(SDL_INIT_VIDEO) != 0) {
        fprintf(stderr, "SDL: %s\n", SDL_GetError());
        return 1;
    }
    bool headless = shot_path || max_frames > 0;
    Uint32 wflags = SDL_WINDOW_ALLOW_HIGHDPI | (headless ? SDL_WINDOW_HIDDEN : SDL_WINDOW_SHOWN);
    SDL_Window *win = SDL_CreateWindow(
        "Parabolic Vessel \xe2\x80\x94 Hard Sphere Gas  [v2.1 cell-list MD]",
        SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, CFG_WIN_W, cfg.sim_h, wflags);
    SDL_Renderer *rend = SDL_CreateRenderer(win, -1,
        headless ? SDL_RENDERER_SOFTWARE
                 : (SDL_RENDERER_ACCELERATED | SDL_RENDERER_PRESENTVSYNC));
    if (!win || !rend) {
        fprintf(stderr, "SDL: %s\n", SDL_GetError());
        return 1;
    }
    // HiDPI/Retina: draw in logical pixels, SDL scales to the backing store
    SDL_RenderSetLogicalSize(rend, CFG_WIN_W, cfg.sim_h);
    g_sprite = make_sprite(rend);
    if (!g_sprite) { fprintf(stderr, "SDL: %s\n", SDL_GetError()); return 1; }

    int n = cfg.n_particles;
    init_particles(n, true);
    sliders_init();
    const double FIXED_DT = 1.0 / 120.0;

    if (shot_path) {
        for (int s = 0; s < shot_steps; s++) {
            step_sim(n, FIXED_DT);
            trail_push(n);
            if (s % 10 == 9) sample_diagnostics(n);
        }
        render_frame(rend, n);
        int ow, oh; SDL_GetRendererOutputSize(rend, &ow, &oh);
        SDL_Surface *sf = SDL_CreateRGBSurfaceWithFormat(0, ow, oh, 24, SDL_PIXELFORMAT_RGB24);
        SDL_RenderReadPixels(rend, NULL, SDL_PIXELFORMAT_RGB24, sf->pixels, sf->pitch);
        SDL_SaveBMP(sf, shot_path);
        SDL_FreeSurface(sf);
        fprintf(stderr, "Screenshot: %s (%dx%d)\n", shot_path, ow, oh);
        SDL_DestroyRenderer(rend); SDL_DestroyWindow(win); SDL_Quit();
        return 0;
    }

    Uint32 t_prev = SDL_GetTicks(), t_last_samp = t_prev;
    double accumulator = 0.0;
    bool running = true;
    long frames = 0;

    while (running) {
        if (max_frames > 0 && ++frames > max_frames) break;
        SDL_Event ev;
        while (SDL_PollEvent(&ev)) {
            if (ev.type == SDL_QUIT) running = false;

            if (ev.type == SDL_MOUSEBUTTONDOWN && ev.button.button == SDL_BUTTON_LEFT) {
                int hit = slider_hit(ev.button.x, ev.button.y, g_slider_panel_top);
                if (hit >= 0) {
                    active_slider = hit;
                    if (slider_set(hit, ev.button.x, g_slider_panel_top)) {
                        n = cfg.n_particles; init_particles(n, true);
                    }
                }
            }
            if (ev.type == SDL_MOUSEBUTTONUP && ev.button.button == SDL_BUTTON_LEFT)
                active_slider = -1;
            if (ev.type == SDL_MOUSEMOTION && active_slider >= 0) {
                if (slider_set(active_slider, ev.motion.x, g_slider_panel_top)) {
                    n = cfg.n_particles; init_particles(n, true);
                }
            }

            if (ev.type == SDL_KEYDOWN)
                switch (ev.key.keysym.sym) {
                case SDLK_ESCAPE: case SDLK_q: running = false; break;
                case SDLK_r: init_particles(n, false); break;
                case SDLK_m: init_particles(n, true);  break;
                case SDLK_t: msd_reset(n); trail_head = trail_count = 0; break;
                case SDLK_c: color_mode = (color_mode + 1) % N_COLOR_MODES; break;
                case SDLK_EQUALS: case SDLK_PLUS:
                    if (n < MAX_PARTICLES) { spawn_particle(n); n++; cfg.n_particles = n; }
                    break;
                case SDLK_MINUS:
                    if (n > 2) { n--; cfg.n_particles = n; if (tracer_i >= n) tracer_i = 0; }
                    break;
                case SDLK_v:
                    if (!rec_active) { record_start(rend); rec_blink = SDL_GetTicks(); }
                    else record_stop();
                    break;
                }
        }

        Uint32 t_now = SDL_GetTicks();
        double frame_time = (t_now - t_prev) / 1000.0;
        if (frame_time > 0.05) frame_time = 0.05;
        t_prev = t_now;
        accumulator += frame_time;
        if (frame_time > 0)
            perf_fps = (perf_fps == 0) ? 1.0 / frame_time : 0.95 * perf_fps + 0.05 / frame_time;

        Uint64 pc0 = SDL_GetPerformanceCounter();
        while (accumulator >= FIXED_DT) {
            step_sim(n, FIXED_DT);
            accumulator -= FIXED_DT;
        }
        double ms = 1000.0 * (double)(SDL_GetPerformanceCounter() - pc0) / SDL_GetPerformanceFrequency();
        perf_step_ms = 0.9 * perf_step_ms + 0.1 * ms;
        trail_push(n);

        if (t_now - t_last_samp >= SAMPLE_MS) {
            sample_diagnostics(n);
            t_last_samp = t_now;
        }

        render_frame(rend, n);
        SDL_RenderPresent(rend);
        capture_frame(rend);
    }

    record_stop();
    SDL_DestroyTexture(g_sprite);
    SDL_DestroyRenderer(rend);
    SDL_DestroyWindow(win);
    SDL_Quit();
    return 0;
}
