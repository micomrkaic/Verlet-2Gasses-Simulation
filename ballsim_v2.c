// ballsim_v2.c — Parabolic Vessel · Hard Sphere Gas  (v2.0)
//
// Physics engine: fixed-dt Verlet integration + cell-list O(N) collision
// detection, replacing the O(N²) event-driven loop from v1.  Supports up to
// 5 000 particles and uses OpenMP to parallelise the Verlet step.
//
// Algorithm per outer fixed-dt step:
//   for each sub-step (N_SUBSTEPS per outer dt):
//     1. Verlet integrate under gravity                [O(N), OpenMP parallel]
//     2. Enforce boundary conditions (walls/ceiling/parabola)  [O(N)]
//     3. Build cell list                               [O(N)]
//     4. Detect & resolve particle-particle overlaps   [O(N·k), k≈9 neighbours]
//        repeated N_RESOLVE_PASSES times for dense packings
//
// Build:
//   gcc -O3 -std=c11 -fopenmp -o ballsim_v2 ballsim_v2.c
//       $(pkg-config --cflags --libs sdl2) -lm
//
// v1 (event-driven, up to 500 particles) is kept intact in
//   verlet_collisions_two_gasses_video.c.

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
#define MAX_PARTICLES     5000  // v2: 10× the v1 cap
#define ENERGY_HIST        280  // samples in energy ring-buffer
#define HIST_BINS           30  // bins in KE-distribution histogram
#define PANEL_W            280  // stats panel width (fixed)
#define N_SUBSTEPS           8  // Verlet sub-steps per outer fixed-dt call
#define N_RESOLVE_PASSES     4  // overlap-resolution passes per sub-step
#define CELL_GRID_MAX      512  // max cells per axis; 512² × 4 B = 1 MB

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
    .radius      = 4.0,        // slightly smaller default to fit more particles
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

static double CFG_PARA_BASE;
static int    CFG_WIN_W;

static void config_derive(void) {
    CFG_PARA_BASE = (cfg.para_a > 1e-9) ? cfg.sim_h - 80.0 : cfg.sim_h;
    CFG_WIN_W     = cfg.sim_w + PANEL_W;
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
        else fprintf(stderr, "Config: unknown key '%s'\n", key);
    }
    fclose(f);
    // Clamp
    if (cfg.n_particles < 1)          cfg.n_particles = 1;
    if (cfg.n_particles > MAX_PARTICLES) cfg.n_particles = MAX_PARTICLES;
    if (cfg.gravity < 0)              cfg.gravity = 0;
    if (cfg.radius  < 1)              cfg.radius  = 1;
    if (cfg.radius  > 20)             cfg.radius  = 20;
    if (cfg.restitution < 0)          cfg.restitution = 0;
    if (cfg.restitution > 1)          cfg.restitution = 1;
    if (cfg.mono_speed  < 1)          cfg.mono_speed  = 1;
    if (cfg.para_a < 0)               cfg.para_a = 0;
    if (cfg.para_a > 0.01)            cfg.para_a = 0.01;
    if (cfg.sim_w  < 200)             cfg.sim_w  = 200;
    if (cfg.sim_w  > 2560)            cfg.sim_w  = 2560;
    if (cfg.sim_h  < 200)             cfg.sim_h  = 200;
    if (cfg.sim_h  > 1440)            cfg.sim_h  = 1440;
    if (cfg.ceil_amp  < 0)            cfg.ceil_amp  = 0;
    if (cfg.ceil_amp  > 200)          cfg.ceil_amp  = 200;
    if (cfg.ceil_freq < 0)            cfg.ceil_freq = 0;
    if (cfg.heavy_frac < 0)           cfg.heavy_frac = 0;
    if (cfg.heavy_frac > 1)           cfg.heavy_frac = 1;
    if (cfg.mass_ratio < 1)           cfg.mass_ratio = 1;
    if (cfg.mass_ratio > 20)          cfg.mass_ratio = 20;
}

static void write_default_config(const char *path) {
    FILE *f = fopen(path, "w");
    if (!f) { fprintf(stderr, "Config: cannot write '%s'\n", path); return; }
    fprintf(f,
        "# Parabolic Vessel v2.0 — config file\n"
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
        "mass_ratio  = %.1f\n",
        cfg.sim_w, cfg.sim_h, cfg.n_particles, cfg.gravity,
        cfg.radius, cfg.restitution, cfg.mono_speed, cfg.para_a,
        cfg.ceil_amp, cfg.ceil_freq, cfg.heavy_frac, cfg.mass_ratio);
    fclose(f);
    fprintf(stderr, "Config: wrote defaults to '%s'\n", path);
}

// ── Parabola vessel ───────────────────────────────────────────────────────────
static inline double para_y(double x) {
    double d = x - cfg.sim_w * 0.5;
    return CFG_PARA_BASE - cfg.para_a * d * d;
}
static inline double para_dydx(double x) {
    return -2.0 * cfg.para_a * (x - cfg.sim_w * 0.5);
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

// ── Video recording (ffmpeg pipe) ─────────────────────────────────────────────
#define REC_FPS 60
static FILE  *rec_pipe    = NULL;
static bool   rec_active  = false;
static Uint32 rec_blink   = 0;
static long   rec_frames  = 0;
static Uint8 *rec_pixels  = NULL;
static char   rec_filename[64];
static Uint32 rec_wall_next = 0;

static bool record_start(int win_w, int win_h) {
    if (rec_active) return true;
    time_t now = time(NULL);
    strftime(rec_filename, sizeof(rec_filename),
             "ballsim_v2_%Y%m%d_%H%M%S.mp4", localtime(&now));
    char cmd[512];
    snprintf(cmd, sizeof(cmd),
        "ffmpeg -y -f rawvideo -pixel_format rgb24"
        " -video_size %dx%d -framerate %d -i pipe:0"
        " -c:v libx264 -preset fast -crf 18 -pix_fmt yuv420p"
        " \"%s\" 2>/dev/null",
        win_w, win_h, REC_FPS, rec_filename);
    rec_pipe = popen(cmd, "w");
    if (!rec_pipe) { fprintf(stderr, "Record: ffmpeg not found?\n"); return false; }
    free(rec_pixels);
    rec_pixels = (Uint8 *)malloc((size_t)win_w * win_h * 3);
    if (!rec_pixels) { pclose(rec_pipe); rec_pipe = NULL; return false; }
    rec_active = true; rec_frames = 0; rec_wall_next = SDL_GetTicks();
    fprintf(stderr, "Record: started → %s\n", rec_filename);
    return true;
}
static void record_stop(void) {
    if (!rec_active) return;
    rec_active = false;
    if (rec_pipe) { pclose(rec_pipe); rec_pipe = NULL; }
    fprintf(stderr, "Record: %ld frames → %s\n", rec_frames, rec_filename);
}
static void capture_frame(SDL_Renderer *rend, int win_w, int win_h) {
    if (!rec_active || !rec_pipe || !rec_pixels) return;
    Uint32 now = SDL_GetTicks();
    if (now < rec_wall_next) return;
    SDL_RenderReadPixels(rend, NULL, SDL_PIXELFORMAT_RGB24, rec_pixels, win_w * 3);
    fwrite(rec_pixels, 1, (size_t)win_w * win_h * 3, rec_pipe);
    rec_frames++;
    rec_wall_next += (Uint32)(1000 / REC_FPS);
    if (now > rec_wall_next + (Uint32)(1000 / REC_FPS)) rec_wall_next = now;
}

// ── Particle ──────────────────────────────────────────────────────────────────
typedef struct {
    double x, y, vx, vy;
    double mass;
    double radius;
    int    species; // 0 = light, 1 = heavy
    Uint8  r, g, b; // species colors (used in SPECIES mode)
    double hue;     // [0,360) — used in RAINBOW mode
} Particle;

static Particle p[MAX_PARTICLES];

// ── Diagnostics buffers ───────────────────────────────────────────────────────
static double eng_buf[ENERGY_HIST];
static int    eng_head = 0, eng_count = 0;
static double eng_max = 1.0, eng_min = 0.0, eng_e0 = 0.0;
static double hist_smooth[HIST_BINS];
static double hist_mb[HIST_BINS];
static double hist_ke_max = 1.0;
static long   collision_count = 0;
static double species_ke[2] = {0, 0};
static int    species_n[2]  = {0, 0};

static double randf(void) { return (double)rand() / ((double)RAND_MAX + 1.0); }

// ── Color cycling ─────────────────────────────────────────────────────────────
#define N_COLOR_MODES 5
static int color_mode = 0;
static const char *color_mode_names[N_COLOR_MODES] = {
    "SPECIES", "SPEED", "KE", "HEIGHT", "RAINBOW"
};

// HSV → RGB  (h in [0,360), s and v in [0,1])
static void hsv_to_rgb(double h, double s, double v,
                        Uint8 *r, Uint8 *g, Uint8 *b) {
    if (s < 1e-9) { *r = *g = *b = (Uint8)(v*255); return; }
    int    i  = (int)(h / 60.0) % 6;
    double f  = h / 60.0 - floor(h / 60.0);
    double p  = v * (1.0 - s);
    double q  = v * (1.0 - f * s);
    double t  = v * (1.0 - (1.0 - f) * s);
    double rv, gv, bv;
    switch (i) {
        case 0: rv=v; gv=t; bv=p; break;
        case 1: rv=q; gv=v; bv=p; break;
        case 2: rv=p; gv=v; bv=t; break;
        case 3: rv=p; gv=q; bv=v; break;
        case 4: rv=t; gv=p; bv=v; break;
        default: rv=v; gv=p; bv=q; break;
    }
    *r = (Uint8)(rv*255); *g = (Uint8)(gv*255); *b = (Uint8)(bv*255);
}

// Compute display color for particle i given current color_mode.
// max_val: frame-level max speed (SPEED) or max KE (KE), used for normalisation.
static void particle_rgb(int i, double max_val,
                          Uint8 *r, Uint8 *g, Uint8 *b) {
    switch (color_mode) {
    case 0: // SPECIES — stored at init
        *r = p[i].r; *g = p[i].g; *b = p[i].b;
        break;
    case 1: { // SPEED — blue=slow, red=fast
        double speed = sqrt(p[i].vx*p[i].vx + p[i].vy*p[i].vy);
        double t     = (max_val > 0) ? speed / max_val : 0.0;
        if (t > 1) t = 1;
        hsv_to_rgb(240.0 * (1.0 - t), 1.0, 1.0, r, g, b);
        break;
    }
    case 2: { // KE — blue=cold, red=hot
        double ke = 0.5 * p[i].mass * (p[i].vx*p[i].vx + p[i].vy*p[i].vy);
        double t  = (max_val > 0) ? ke / max_val : 0.0;
        if (t > 1) t = 1;
        hsv_to_rgb(240.0 * (1.0 - t), 1.0, 1.0, r, g, b);
        break;
    }
    case 3: { // HEIGHT — rainbow by y position (red=top, blue=bottom)
        double t = p[i].y / cfg.sim_h;
        if (t < 0) t = 0;
        if (t > 1) t = 1;
        hsv_to_rgb(t * 300.0, 1.0, 1.0, r, g, b);
        break;
    }
    case 4: // RAINBOW — fixed hue assigned at init
        hsv_to_rgb(p[i].hue, 1.0, 1.0, r, g, b);
        break;
    default:
        *r = p[i].r; *g = p[i].g; *b = p[i].b;
    }
}

// ── Maxwell-Boltzmann curve (2D: f(KE) = exp(-KE/kT)) ───────────────────────
static void compute_mb_curve(double mean_ke) {
    if (mean_ke < 1.0) mean_ke = 1.0;
    double kT = mean_ke, peak = 0.0;
    double raw[HIST_BINS];
    for (int b = 0; b < HIST_BINS; b++) {
        double ke = (b + 0.5) / HIST_BINS * hist_ke_max;
        raw[b] = exp(-ke / kT);
        if (raw[b] > peak) peak = raw[b];
    }
    for (int b = 0; b < HIST_BINS; b++)
        hist_mb[b] = (peak > 0) ? raw[b] / peak : 0.0;
}

// ── Particle initialisation ───────────────────────────────────────────────────
static void init_particles(int n, bool monoenergetic) {
    double max_r = cfg.radius;
    for (int i = 0; i < n; i++) {
        int is_heavy  = (randf() < cfg.heavy_frac) ? 1 : 0;
        p[i].species  = is_heavy;
        p[i].mass     = is_heavy ? cfg.mass_ratio : 1.0;
        p[i].radius   = cfg.radius * cbrt(p[i].mass);
        if (p[i].radius > max_r) max_r = p[i].radius;

        double angle = randf() * 2.0 * M_PI;
        double speed = monoenergetic ? cfg.mono_speed : 80.0 + randf() * 300.0;
        double v     = speed / sqrt(p[i].mass);
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
        p[i].hue = randf() * 360.0;
    }
    // Grid placement
    double gap    = 2.0;
    double cell   = 2.0 * max_r + gap;
    double piston = cfg.radius + cfg.ceil_amp;
    double y_start= piston + max_r + gap;
    double x_lo   = max_r + gap;
    double x_hi   = cfg.sim_w - max_r - gap;
    int cols = (int)((x_hi - x_lo) / cell);
    if (cols < 1) cols = 1;
    for (int i = 0; i < n; i++) {
        int col = i % cols, row = i / cols;
        p[i].x = x_lo + col * cell + max_r + (randf() - 0.5) * gap;
        p[i].y = y_start + row * cell             + (randf() - 0.5) * gap;
        double ri = p[i].radius;
        if (p[i].x < ri + 1)             p[i].x = ri + 1;
        if (p[i].x > cfg.sim_w - ri - 1) p[i].x = cfg.sim_w - ri - 1;
        if (p[i].y < ri + 1)             p[i].y = ri + 1;
        if (p[i].y > cfg.sim_h - ri - 1) p[i].y = cfg.sim_h - ri - 1;
    }
    eng_head = eng_count = 0;
    eng_max = 1.0; eng_min = 0.0; eng_e0 = 0.0;
    for (int b = 0; b < HIST_BINS; b++) hist_smooth[b] = hist_mb[b] = 0.0;
    hist_ke_max = 1.0;
    collision_count = 0;
    sim_time = 0.0;
}

// ── Energy ────────────────────────────────────────────────────────────────────
static double compute_energy(int n, double *ke_out, double *pe_out) {
    double ke = 0, pe = 0;
    species_ke[0] = species_ke[1] = 0;
    species_n[0]  = species_n[1]  = 0;
    for (int i = 0; i < n; i++) {
        double pke = 0.5 * p[i].mass * (p[i].vx*p[i].vx + p[i].vy*p[i].vy);
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
    double norm = (eng_e0 != 0.0) ? e / eng_e0 : 1.0;
    eng_buf[eng_head] = norm;
    eng_head = (eng_head + 1) % ENERGY_HIST;
    if (eng_count < ENERGY_HIST) eng_count++;
    eng_max = eng_min = eng_buf[0];
    for (int i = 0; i < eng_count; i++) {
        if (eng_buf[i] > eng_max) eng_max = eng_buf[i];
        if (eng_buf[i] < eng_min) eng_min = eng_buf[i];
    }
    if (eng_max < 1.0) eng_max = 1.0;
    if (eng_min > 1.0) eng_min = 1.0;
    if (eng_max - eng_min < 0.002) { eng_max = 1.001; eng_min = 0.999; }
}

static void update_histogram(int n) {
    double max_ke = 1.0;
    for (int i = 0; i < n; i++) {
        double ke = 0.5*p[i].mass*(p[i].vx*p[i].vx + p[i].vy*p[i].vy);
        if (ke > max_ke) max_ke = ke;
    }
    hist_ke_max += (max_ke * 1.1 - hist_ke_max) * 0.05;
    double raw[HIST_BINS] = {0};
    for (int i = 0; i < n; i++) {
        double ke = 0.5*p[i].mass*(p[i].vx*p[i].vx + p[i].vy*p[i].vy);
        int b = (int)(ke / hist_ke_max * HIST_BINS);
        if (b >= HIST_BINS) b = HIST_BINS - 1;
        raw[b] += 1.0;
    }
    double peak = 0;
    for (int b = 0; b < HIST_BINS; b++) if (raw[b] > peak) peak = raw[b];
    if (peak > 0)
        for (int b = 0; b < HIST_BINS; b++) raw[b] /= peak;
    for (int b = 0; b < HIST_BINS; b++)
        hist_smooth[b] = hist_smooth[b] * 0.85 + raw[b] * 0.15;
    double ke_total = 0;
    for (int i = 0; i < n; i++)
        ke_total += 0.5*p[i].mass*(p[i].vx*p[i].vx + p[i].vy*p[i].vy);
    compute_mb_curve(n > 0 ? ke_total / n : 1.0);
}

// ═══════════════════════════════════════════════════════════════════════════════
// ── v2 Physics: cell-list Verlet MD ──────────────────────────────────────────
// ═══════════════════════════════════════════════════════════════════════════════

// Cell list: singly-linked lists, one per grid cell.
// cell_head[cidx] = first particle index in cell cidx, or -1 if empty.
// cell_next[i]    = next particle in the same cell as p[i], or -1.
// dirty[]         = indices of cells that were non-empty last build (so we can
//                   clear only those, not the whole 1 MB array, each sub-step).
static int g_cell_head[CELL_GRID_MAX * CELL_GRID_MAX];
static int g_cell_next[MAX_PARTICLES];
static int g_dirty[MAX_PARTICLES];
static int g_ndirty;
static int    g_ncx, g_ncy; // current grid dimensions
static double g_inv_cs;     // 1 / cell_size

// Call once before the first simulation step.
static void cell_list_init(void) {
    memset(g_cell_head, 0xFF, sizeof(g_cell_head)); // 0xFF bytes → -1 as int32
    g_ndirty = 0;
}

static void build_cells(int n) {
    // Largest possible particle radius determines cell size.
    // Cell size ≥ 2*max_r ensures that any two overlapping particles are in the
    // same or adjacent cells, so the 3×3 neighbourhood check never misses a pair.
    double max_r = cfg.radius * cbrt(cfg.mass_ratio > 1.0 ? cfg.mass_ratio : 1.0);
    double cs    = 2.0 * max_r + 0.5; // +0.5 px safety margin
    g_inv_cs = 1.0 / cs;
    g_ncx = (int)(cfg.sim_w * g_inv_cs) + 2;
    g_ncy = (int)(cfg.sim_h * g_inv_cs) + 2;
    if (g_ncx > CELL_GRID_MAX) g_ncx = CELL_GRID_MAX;
    if (g_ncy > CELL_GRID_MAX) g_ncy = CELL_GRID_MAX;

    // Clear only cells touched in the previous build — O(N) not O(grid).
    for (int k = 0; k < g_ndirty; k++)
        g_cell_head[g_dirty[k]] = -1;
    g_ndirty = 0;

    for (int i = 0; i < n; i++) {
        int cx = (int)(p[i].x * g_inv_cs);
        int cy = (int)(p[i].y * g_inv_cs);
        if (cx < 0) cx = 0;
        if (cx >= g_ncx) cx = g_ncx - 1;
        if (cy < 0) cy = 0;
        if (cy >= g_ncy) cy = g_ncy - 1;
        int cidx = cy * g_ncx + cx;
        if (g_cell_head[cidx] == -1)        // first particle in this cell
            g_dirty[g_ndirty++] = cidx;     // mark for cleanup next build
        g_cell_next[i]   = g_cell_head[cidx];
        g_cell_head[cidx] = i;
    }
}

// Boundary conditions: walls, ceiling (piston), parabola/flat floor.
// Called with the absolute time at the END of this sub-step.
static void apply_boundaries(int n, double t_abs) {
    double cy  = ceil_y(t_abs);
    double cvy = ceil_vy(t_abs);

    for (int i = 0; i < n; i++) {
        double ri = p[i].radius;

        // ── Ceiling (piston) ─────────────────────────────────────────────────
        if (p[i].y - ri < cy) {
            p[i].y = cy + ri + 1e-4;
            double rel = p[i].vy - cvy;
            if (rel < 0) p[i].vy = -rel * cfg.restitution + cvy;
        }

        // ── Side walls (only for particles above the parabola rim) ───────────
        if (cfg.para_a < 1e-9 || p[i].y < para_y(ri)) {
            if (p[i].x - ri < 0) {
                p[i].x = ri + 1e-4;
                if (p[i].vx < 0) p[i].vx = -p[i].vx * cfg.restitution;
            }
            if (p[i].x + ri > cfg.sim_w) {
                p[i].x = cfg.sim_w - ri - 1e-4;
                if (p[i].vx > 0) p[i].vx = -p[i].vx * cfg.restitution;
            }
        }

        // ── Floor ────────────────────────────────────────────────────────────
        if (cfg.para_a > 1e-9) {
            // Parabola: reflect about surface normal.
            double floor = para_y(p[i].x);
            if (p[i].y + ri > floor) {
                p[i].y = floor - ri - 1e-4;
                // Outward normal (points INTO floor, i.e. downward at centre):
                // tangent = (1, dydx), outward normal = (-dydx, 1)/|n|
                double dydx = para_dydx(p[i].x);
                double nx   = -dydx, ny = 1.0;
                double nlen = sqrt(nx*nx + ny*ny);
                nx /= nlen; ny /= nlen;
                double vn = p[i].vx * nx + p[i].vy * ny;
                if (vn > 0) { // moving into floor
                    p[i].vx = (p[i].vx - 2.0 * vn * nx) * cfg.restitution;
                    p[i].vy = (p[i].vy - 2.0 * vn * ny) * cfg.restitution;
                }
            }
        } else {
            // Flat floor
            double floor = CFG_PARA_BASE - ri;
            if (p[i].y > floor) {
                p[i].y = floor - 1e-4;
                if (p[i].vy > 0) p[i].vy = -p[i].vy * cfg.restitution;
            }
        }
    }
}

// One pass of particle-particle overlap detection + impulse resolution.
// Each pair (i,j) with i<j is processed at most once per pass.
// Uses the cell list built by build_cells() — do not modify positions outside
// a 2*cell_size neighbourhood or pairs may be missed.
static void resolve_collisions(int n) {
    for (int i = 0; i < n; i++) {
        int cx = (int)(p[i].x * g_inv_cs);
        int cy = (int)(p[i].y * g_inv_cs);
        if (cx < 0) cx = 0;
        if (cx >= g_ncx) cx = g_ncx - 1;
        if (cy < 0) cy = 0;
        if (cy >= g_ncy) cy = g_ncy - 1;

        for (int dcy = -1; dcy <= 1; dcy++) {
            int ncy2 = cy + dcy;
            if (ncy2 < 0 || ncy2 >= g_ncy) continue;
            for (int dcx = -1; dcx <= 1; dcx++) {
                int ncx2 = cx + dcx;
                if (ncx2 < 0 || ncx2 >= g_ncx) continue;
                int cidx = ncy2 * g_ncx + ncx2;

                for (int j = g_cell_head[cidx]; j >= 0; j = g_cell_next[j]) {
                    if (j <= i) continue; // process each pair exactly once

                    double dx  = p[j].x - p[i].x;
                    double dy  = p[j].y - p[i].y;
                    double d2  = dx*dx + dy*dy;
                    double tch = p[i].radius + p[j].radius;
                    if (d2 >= tch * tch || d2 < 1e-20) continue;

                    double d  = sqrt(d2);
                    double nx = dx / d, ny = dy / d;

                    // ── Velocity impulse (elastic, with restitution) ─────────
                    double dvn = (p[j].vx - p[i].vx)*nx + (p[j].vy - p[i].vy)*ny;
                    if (dvn < 0) { // particles approaching
                        double mi = p[i].mass, mj = p[j].mass;
                        double J  = -(1.0 + cfg.restitution) * mi * mj / (mi + mj) * dvn;
                        p[i].vx -= (J / mi) * nx;
                        p[i].vy -= (J / mi) * ny;
                        p[j].vx += (J / mj) * nx;
                        p[j].vy += (J / mj) * ny;
                        collision_count++;
                    }

                    // ── Position separation: push apart along contact normal ──
                    // Split by mass ratio so the lighter particle moves more.
                    double overlap = tch - d;
                    if (overlap <= 0) continue;
                    double push = overlap * 0.5 + 1e-5;
                    double mi = p[i].mass, mj = p[j].mass;
                    double tot = mi + mj;
                    p[i].x -= push * (mj / tot) * nx;
                    p[i].y -= push * (mj / tot) * ny;
                    p[j].x += push * (mi / tot) * nx;
                    p[j].y += push * (mi / tot) * ny;
                }
            }
        }
    }
}

// ── Main simulation step ──────────────────────────────────────────────────────
// Advances simulation by dt seconds, internally using N_SUBSTEPS sub-steps.
// Called from the outer fixed-dt accumulator loop at 120 Hz.
static void step_sim(int n, double dt) {
    double sub_dt = dt / N_SUBSTEPS;

    for (int sub = 0; sub < N_SUBSTEPS; sub++) {
        double t_end = sim_time + (sub + 1) * sub_dt;

        // ── 1. Verlet integrate under gravity ─────────────────────────────────
        // Embarrassingly parallel: no shared writes.
#ifdef _OPENMP
#       pragma omp parallel for schedule(static)
#endif
        for (int i = 0; i < n; i++) {
            p[i].x  += p[i].vx * sub_dt;
            p[i].y  += p[i].vy * sub_dt + 0.5 * cfg.gravity * sub_dt * sub_dt;
            p[i].vy += cfg.gravity * sub_dt;
        }

        // ── 2. Boundary conditions ────────────────────────────────────────────
        apply_boundaries(n, t_end);

        // ── 3. Cell list ──────────────────────────────────────────────────────
        build_cells(n);

        // ── 4. Collision resolution (N passes for dense packings) ─────────────
        for (int pass = 0; pass < N_RESOLVE_PASSES; pass++) {
            resolve_collisions(n);
            // Re-clamp positions to hard walls after each resolution pass
            // (position corrections can push particles outside boundaries).
            for (int i = 0; i < n; i++) {
                double ri = p[i].radius;
                if (p[i].x < ri)             p[i].x = ri;
                if (p[i].x > cfg.sim_w - ri) p[i].x = cfg.sim_w - ri;
                double cy2 = ceil_y(t_end);
                if (p[i].y - ri < cy2)       p[i].y = cy2 + ri;
                if (cfg.para_a < 1e-9) {
                    if (p[i].y > CFG_PARA_BASE - ri) p[i].y = CFG_PARA_BASE - ri;
                } else {
                    if (p[i].y + ri > para_y(p[i].x))
                        p[i].y = para_y(p[i].x) - ri;
                }
            }
        }
    }
    sim_time += dt;
}

// ═══════════════════════════════════════════════════════════════════════════════
// ── Drawing ───────────────────────────────────────────────────────────────────
// ═══════════════════════════════════════════════════════════════════════════════

// Fast scanline circle — 2R+1 SDL_RenderDrawLine calls instead of ~200 point
// calls, essential for 5000 particles at 60 fps.
static void draw_particle(SDL_Renderer *rend, double x, double y, double rad,
                          Uint8 r, Uint8 g, Uint8 b) {
    int cx = (int)x, cy2 = (int)y;
    int R  = (int)(rad + 0.5);
    if (R < 1) R = 1;
    for (int dy = -R; dy <= R; dy++) {
        int w = (int)sqrt((double)(R*R - dy*dy));
        if (w == 0) w = 1;
        // Simple radial shading: brighter at centre
        double br = 1.0 - 0.38 * (double)(dy*dy) / (double)(R*R + 1);
        SDL_SetRenderDrawColor(rend, (Uint8)(r*br), (Uint8)(g*br), (Uint8)(b*br), 255);
        SDL_RenderDrawLine(rend, cx - w, cy2 + dy, cx + w, cy2 + dy);
    }
    // Small specular highlight (top-left)
    SDL_SetRenderDrawColor(rend,
        (Uint8)fmin(255, r * 1.6),
        (Uint8)fmin(255, g * 1.6),
        (Uint8)fmin(255, b * 1.6), 255);
    SDL_RenderDrawPoint(rend, cx - R/3, cy2 - R/3);
}

static void draw_vessel(SDL_Renderer *rend) {
    SDL_SetRenderDrawColor(rend, 0x55, 0x77, 0xAA, 0xFF);
    for (int xs = 0; xs < cfg.sim_w - 1; xs++) {
        int y0 = (int)para_y(xs), y1 = (int)para_y(xs + 1);
        for (int t = -1; t <= 1; t++)
            SDL_RenderDrawLine(rend, xs, y0+t, xs+1, y1+t);
    }
    {
        int cy2 = (int)ceil_y(sim_time);
        SDL_SetRenderDrawColor(rend, 0x1A, 0x22, 0x3A, 0xFF);
        SDL_Rect piston = {0, 0, cfg.sim_w, cy2};
        SDL_RenderFillRect(rend, &piston);
        SDL_SetRenderDrawColor(rend, 0x55, 0x77, 0xAA, 0xFF);
        for (int t = -1; t <= 1; t++)
            SDL_RenderDrawLine(rend, 0, cy2+t, cfg.sim_w, cy2+t);
    }
    double rim_y = (cfg.para_a > 1e-9) ? para_y(cfg.radius) : CFG_PARA_BASE;
    SDL_Rect lw = {0, 0, 2, (int)rim_y}, rw = {cfg.sim_w-2, 0, 2, (int)rim_y};
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
    {0x0,0x0,0x0,0x0,0x2},{0x0,0x0,0x7,0x0,0x0},{0x1,0x1,0x2,0x4,0x4},
    {0x3,0x4,0x4,0x4,0x3},{0x6,0x1,0x1,0x1,0x6},{0x0,0x0,0x0,0x0,0x0},
};
static int gl(char c) {
    if (c>='0'&&c<='9') return c-'0';
    if (c>='A'&&c<='Z') return 10+(c-'A');
    if (c>='a'&&c<='z') return 10+(c-'a');
    if (c=='.') return 36;
    if (c=='-') return 37;
    if (c=='/') return 38;
    if (c=='(') return 39;
    if (c==')') return 40;
    return 41;
}
static void ds(SDL_Renderer *rend, const char *s, int x, int y, int sc,
               Uint8 r, Uint8 g, Uint8 b) {
    SDL_SetRenderDrawColor(rend, r, g, b, 255);
    for (int ci = 0; s[ci]; ci++)
        for (int row = 0; row < 5; row++)
            for (int col = 0; col < 3; col++)
                if (F[gl(s[ci])][row] & (0x4 >> col)) {
                    SDL_Rect px = {x + ci*(3*sc+sc) + col*sc, y + row*sc, sc, sc};
                    SDL_RenderFillRect(rend, &px);
                }
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
    SI("PARTICLES",  1,  5000, n_particles, 1); // v2: 10× range
    SD("GRAVITY",    0,  2000, gravity,     0);
    SD("RADIUS",     1,    20, radius,      1);
    SD("RESTITUTION",0,     1, restitution, 0);
    SD("MONO SPEED", 10, 1000, mono_speed,  1);
    SD("PARA A",     0,  0.01, para_a,      1);
    SD("PISTON AMP", 0,   150, ceil_amp,    0);
    SD("PISTON FREQ",0,    10, ceil_freq,   0);
    SD("HEAVY FRAC", 0,     1, heavy_frac,  1);
    SD("MASS RATIO", 1,    20, mass_ratio,  1);
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
    return (SDL_Rect){cfg.sim_w+8, top+18+i*30+16, PANEL_W-16, 5};
}
static void draw_sliders(SDL_Renderer *rend, int top) {
    ds(rend, "PARAMETERS", cfg.sim_w+8, top, 2, 0xAA, 0xBB, 0xCC);
    for (int i = 0; i < n_sliders; i++) {
        Slider *s = &sliders[i];
        SDL_Rect tr = slider_track(i, top);
        int ly = tr.y - 13;
        ds(rend, s->label, tr.x, ly, 2, 0x88, 0x99, 0xBB);
        char buf[16];
        double v = slider_get(i);
        if (s->ivalue)       snprintf(buf,sizeof(buf),"%d",(int)v);
        else if (s->max<=0.1)snprintf(buf,sizeof(buf),"%.4f",v);
        else if (s->max<=1.0)snprintf(buf,sizeof(buf),"%.3f",v);
        else if (s->max<=10.0)snprintf(buf,sizeof(buf),"%.2f",v);
        else                  snprintf(buf,sizeof(buf),"%.1f",v);
        ds(rend, buf, tr.x+tr.w-(int)strlen(buf)*8, ly, 2, 0xFF, 0xDD, 0x88);
        SDL_SetRenderDrawColor(rend, 0x22, 0x2A, 0x3A, 0xFF);
        SDL_RenderFillRect(rend, &tr);
        double frac = (v-s->min)/(s->max-s->min);
        if (frac<0) frac=0;
        if (frac>1) frac=1;
        SDL_Rect fill = tr; fill.w = (int)(frac*tr.w);
        SDL_SetRenderDrawColor(rend,
            i==active_slider?0x88:0x44,
            i==active_slider?0xCC:0x88,
            i==active_slider?0xFF:0xCC, 0xFF);
        if (fill.w>0) SDL_RenderFillRect(rend, &fill);
        SDL_Rect thumb = {tr.x+fill.w-4, tr.y-4, 8, 13};
        SDL_SetRenderDrawColor(rend, 0xDD, 0xEE, 0xFF, 0xFF);
        SDL_RenderFillRect(rend, &thumb);
    }
}
static int slider_hit(int mx, int my, int top) {
    for (int i = 0; i < n_sliders; i++) {
        SDL_Rect tr = slider_track(i,top);
        SDL_Rect h  = {tr.x, tr.y-6, tr.w, tr.h+12};
        if (mx>=h.x&&mx<h.x+h.w&&my>=h.y&&my<h.y+h.h) return i;
    }
    return -1;
}
static int slider_set(int idx, int mx, int top) {
    Slider *s = &sliders[idx];
    SDL_Rect tr = slider_track(idx, top);
    double frac = (double)(mx-tr.x)/tr.w;
    if (frac<0) frac=0;
        if (frac>1) frac=1;
    slider_put(idx, s->min + frac*(s->max-s->min));
    config_derive();
    return s->needs_reset;
}

// ── Stats row helper ──────────────────────────────────────────────────────────
static int draw_stat(SDL_Renderer *rend, int px, int y,
                     const char *label, const char *val,
                     Uint8 vr, Uint8 vg, Uint8 vb) {
    ds(rend, label, px, y,    2, 0xAA, 0xBB, 0xCC);
    ds(rend, val,   px, y+12, 3, vr,   vg,   vb);
    return y + 30;
}

// ── Stats panel ───────────────────────────────────────────────────────────────
static void draw_panel(SDL_Renderer *rend, int n, double total_e,
                        double ke, double pe) {
    SDL_SetRenderDrawColor(rend, 0x0C, 0x0F, 0x1C, 0xFF);
    SDL_Rect bg = {cfg.sim_w, 0, PANEL_W, cfg.sim_h};
    SDL_RenderFillRect(rend, &bg);
    SDL_SetRenderDrawColor(rend, 0x44, 0x55, 0x77, 0xFF);
    SDL_RenderDrawLine(rend, cfg.sim_w, 0, cfg.sim_w, cfg.sim_h);

    int px = cfg.sim_w + 8, ty = 10;
    char buf[32];

    // REC indicator
    {
        int rx = cfg.sim_w + PANEL_W - 8;
        if (rec_active) {
            bool dot_on = ((SDL_GetTicks()-rec_blink)%900) < 600;
            if (dot_on) {
                SDL_SetRenderDrawColor(rend, 0xEE, 0x22, 0x22, 0xFF);
                for (int dy=-5;dy<=5;dy++) for (int dx=-5;dx<=5;dx++)
                    if (dx*dx+dy*dy<=25) SDL_RenderDrawPoint(rend,rx-6+dx,ty+5+dy);
            }
            ds(rend,"REC",rx-38,ty,2,0xEE,0x33,0x33);
            snprintf(buf,sizeof(buf),"%ldf",rec_frames);
            ds(rend,buf,rx-(int)strlen(buf)*8,ty+14,2,0xAA,0x44,0x44);
        } else {
            SDL_SetRenderDrawColor(rend,0x33,0x33,0x33,0xFF);
            for (int dy=-4;dy<=4;dy++) for (int dx=-4;dx<=4;dx++)
                if (dx*dx+dy*dy<=16) SDL_RenderDrawPoint(rend,rx-6+dx,ty+4+dy);
            ds(rend,"V-REC",rx-46,ty,2,0x44,0x44,0x55);
        }
    }

    // Version badge
    ds(rend, "V2.0", px, ty, 2, 0x33, 0x88, 0xFF);
    ty += 18;

    snprintf(buf,sizeof(buf),"%dx%d",cfg.sim_w,cfg.sim_h);
    ty = draw_stat(rend,px,ty,"BOX SIZE (cfg sim-w/h)",buf,0x88,0x99,0xAA);
    ty += 4;

    snprintf(buf,sizeof(buf),"%d",n);
    ty = draw_stat(rend,px,ty,"PARTICLES",buf,0xFF,0xFF,0xFF);
    ty += 4;

    ty = draw_stat(rend,px,ty,"COLOR MODE  (C)",
                   color_mode_names[color_mode],0xFF,0xDD,0xFF);
    ty += 4;

    snprintf(buf,sizeof(buf),"%.4f",(eng_e0>0)?total_e/eng_e0:1.0);
    ty = draw_stat(rend,px,ty,"TOTAL ENERGY (E/E0)",buf,0xFF,0xCC,0x44);
    ty += 4;

    snprintf(buf,sizeof(buf),"%.4f",(total_e!=0)?ke/total_e:0.0);
    ty = draw_stat(rend,px,ty,"KINETIC (KE/E)",buf,0x44,0xFF,0xAA);
    ty += 4;

    snprintf(buf,sizeof(buf),"%.4f",(total_e!=0)?pe/total_e:0.0);
    ty = draw_stat(rend,px,ty,"POTENTIAL (PE/E)",buf,0x88,0xAA,0xFF);
    ty += 4;

    double mean_ke = (n>0) ? ke/n : 0.0;
    snprintf(buf,sizeof(buf),"%.1f",mean_ke);
    ty = draw_stat(rend,px,ty,"TEMPERATURE (KE/N)",buf,0xFF,0x66,0x44);
    ty += 4;

    {
        double lke = species_n[0]>0 ? species_ke[0]/species_n[0] : 0.0;
        double hke = species_n[1]>0 ? species_ke[1]/species_n[1] : 0.0;
        snprintf(buf,sizeof(buf),"%.1f",lke);
        ty = draw_stat(rend,px,ty,"LIGHT KE/N",buf,0x33,0xCC,0xFF);
        ty += 4;
        snprintf(buf,sizeof(buf),"%.1f",hke);
        ty = draw_stat(rend,px,ty,"HEAVY KE/N",buf,0xFF,0x88,0x33);
        ty += 4;
    }

    snprintf(buf,sizeof(buf),"%ld",collision_count);
    ty = draw_stat(rend,px,ty,"COLLISIONS",buf,0xCC,0x88,0xFF);
    ty += 10;

    // ── Energy trace ──────────────────────────────────────────────────────────
    ds(rend,"ENERGY / E0",px,ty,2,0xAA,0xBB,0xCC);
    ty += 13;
    int gx=cfg.sim_w+4, gw=PANEL_W-8, gh=80, gy=ty;
    SDL_SetRenderDrawColor(rend,0x06,0x08,0x10,0xFF);
    SDL_Rect gb={gx,gy,gw,gh}; SDL_RenderFillRect(rend,&gb);
    double range = eng_max - eng_min;
#define TRACE_Y(v) (gy+gh-1-(int)(((v)-eng_min)/range*(gh-2)))
    SDL_SetRenderDrawBlendMode(rend,SDL_BLENDMODE_BLEND);
    SDL_SetRenderDrawColor(rend,0x22,0x33,0x55,80);
    for (int f=1;f<=3;f++)
        SDL_RenderDrawLine(rend,gx,gy+gh-gh*f/4,gx+gw,gy+gh-gh*f/4);
    int ref_y=TRACE_Y(1.0);
    if (ref_y>=gy&&ref_y<gy+gh) {
        SDL_SetRenderDrawColor(rend,0xFF,0xFF,0xFF,60);
        SDL_RenderDrawLine(rend,gx,ref_y,gx+gw,ref_y);
    }
    SDL_SetRenderDrawBlendMode(rend,SDL_BLENDMODE_NONE);
    if (eng_count>1) {
        int cols=(eng_count<gw)?eng_count:gw;
        for (int col=0;col<cols;col++) {
            int idx=(eng_head-1-col+ENERGY_HIST)%ENERGY_HIST;
            double v=eng_buf[idx];
            int bar_y=TRACE_Y(v);
            if (bar_y<gy) bar_y=gy;
            if (bar_y>=gy+gh) bar_y=gy+gh-1;
            int bx=gx+gw-1-col;
            double dev=fabs(v-1.0)/(range*0.5+1e-9); if (dev>1.0) dev=1.0;
            SDL_SetRenderDrawColor(rend,(Uint8)(dev*255),(Uint8)((1.0-dev)*200),0x18,0xFF);
            SDL_RenderDrawLine(rend,bx,bar_y,bx,gy+gh-1);
        }
    }
    SDL_SetRenderDrawColor(rend,0x2A,0x3A,0x55,0xFF);
    SDL_RenderDrawRect(rend,&gb);
    if (eng_count>0) {
        double cur=eng_buf[(eng_head-1+ENERGY_HIST)%ENERGY_HIST];
        double pct=(cur-1.0)*100.0;
        char devbuf[16]; snprintf(devbuf,sizeof(devbuf),"%+.2f PCT",pct);
        ds(rend,"DRIFT",px,ty+gh+2,2,0xAA,0xBB,0xCC);
        Uint8 dcr=(Uint8)(fabs(pct)>1.0?220:fabs(pct)>0.1?180:100);
        Uint8 dcg=(Uint8)(fabs(pct)>1.0?80:fabs(pct)>0.1?180:220);
        ds(rend,devbuf,px+40,ty+gh+2,2,dcr,dcg,0x44);
    }
#undef TRACE_Y
    ty += gh + 18;

    // ── KE histogram ──────────────────────────────────────────────────────────
    ds(rend,"KE DISTRIBUTION vs MB",px,ty,2,0xAA,0xBB,0xCC);
    ty += 13;
    int hx=cfg.sim_w+4, hw=PANEL_W-8, hh=85, hy=ty;
    SDL_SetRenderDrawColor(rend,0x06,0x08,0x10,0xFF);
    SDL_Rect hb={hx,hy,hw,hh}; SDL_RenderFillRect(rend,&hb);
    SDL_SetRenderDrawBlendMode(rend,SDL_BLENDMODE_BLEND);
    SDL_SetRenderDrawColor(rend,0x22,0x33,0x55,80);
    for (int f=1;f<=3;f++)
        SDL_RenderDrawLine(rend,hx,hy+hh-hh*f/4,hx+hw,hy+hh-hh*f/4);
    SDL_SetRenderDrawBlendMode(rend,SDL_BLENDMODE_NONE);
    double bin_w=(double)hw/HIST_BINS;
    for (int b=0;b<HIST_BINS;b++) {
        int bx=(int)(hx+b*bin_w), bw=(int)bin_w-1; if (bw<1) bw=1;
        int bar_h=(int)(hist_smooth[b]*(hh-2));
        if (bar_h>0) {
            SDL_SetRenderDrawBlendMode(rend,SDL_BLENDMODE_BLEND);
            SDL_SetRenderDrawColor(rend,0x22,0xAA,0xCC,180);
            SDL_Rect br={bx,hy+hh-bar_h,bw,bar_h}; SDL_RenderFillRect(rend,&br);
            SDL_SetRenderDrawBlendMode(rend,SDL_BLENDMODE_NONE);
        }
        int mb_y=(int)(hy+hh-hist_mb[b]*(hh-2));
        if (b>0) {
            int prev_mb_y=(int)(hy+hh-hist_mb[b-1]*(hh-2));
            SDL_SetRenderDrawColor(rend,0xFF,0xDD,0x33,0xFF);
            SDL_RenderDrawLine(rend,
                (int)(hx+(b-0.5)*bin_w),prev_mb_y,
                (int)(hx+(b+0.5)*bin_w),mb_y);
        }
        (void)mb_y;
    }
    {
        int prev_y=(int)(hy+hh-hist_mb[HIST_BINS-2]*(hh-2));
        int last_y=(int)(hy+hh-hist_mb[HIST_BINS-1]*(hh-2));
        SDL_SetRenderDrawColor(rend,0xFF,0xDD,0x33,0xFF);
        SDL_RenderDrawLine(rend,
            (int)(hx+(HIST_BINS-1.5)*bin_w),prev_y,
            (int)(hx+(HIST_BINS-0.5)*bin_w),last_y);
    }
    SDL_SetRenderDrawColor(rend,0x2A,0x3A,0x55,0xFF);
    SDL_RenderDrawRect(rend,&hb);
    ty += hh + 6;

    // Legend
    SDL_SetRenderDrawColor(rend,0x22,0xAA,0xCC,0xFF);
    SDL_Rect lc={px,ty+2,8,8}; SDL_RenderFillRect(rend,&lc);
    ds(rend,"MEASURED",px+11,ty,2,0x44,0xCC,0xEE);
    ty += 14;
    SDL_SetRenderDrawColor(rend,0xFF,0xDD,0x33,0xFF);
    SDL_RenderDrawLine(rend,px,ty+4,px+8,ty+4);
    ds(rend,"MAXWELL-BOLTZMANN",px+11,ty,2,0xDD,0xBB,0x44);
    ty += 20;

    // Species legend
    SDL_SetRenderDrawColor(rend,0x33,0xCC,0xFF,0xFF);
    SDL_Rect sl={px,ty+2,8,8}; SDL_RenderFillRect(rend,&sl);
    ds(rend,"LIGHT",px+11,ty,2,0x33,0xCC,0xFF);
    SDL_SetRenderDrawColor(rend,0xFF,0x88,0x33,0xFF);
    SDL_Rect sh={px+60,ty+2,8,8}; SDL_RenderFillRect(rend,&sh);
    ds(rend,"HEAVY",px+71,ty,2,0xFF,0x88,0x33);
    ty += 14;

    // Divider
    SDL_SetRenderDrawColor(rend,0x2A,0x3A,0x55,0xFF);
    SDL_RenderDrawLine(rend,cfg.sim_w+4,ty,cfg.sim_w+PANEL_W-4,ty);
    ty += 6;

    // Sliders
    g_slider_panel_top = ty;
    draw_sliders(rend, ty);

    // OpenMP thread count
#ifdef _OPENMP
    char omp_buf[24];
    snprintf(omp_buf, sizeof(omp_buf), "OMP THREADS-%d", omp_get_max_threads());
    ds(rend, omp_buf, px, cfg.sim_h-24, 1, 0x33, 0x66, 0x44);
#endif

    ds(rend, "R-RAND  M-MONO  C-COLOR  V-REC  Q-QUIT",
       px, cfg.sim_h-14, 1, 0x55, 0x66, 0x88);
}

// ── Main ──────────────────────────────────────────────────────────────────────
int main(int argc, char **argv) {
    const char *config_path = "ballsim_v2.cfg";
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i],"--config") && i+1<argc) {
            config_path = argv[++i];
        } else if (!strcmp(argv[i],"--write-config")) {
            const char *out = (i+1<argc) ? argv[++i] : "ballsim_v2.cfg";
            config_derive();
            write_default_config(out);
            return 0;
        } else if (!strcmp(argv[i],"--help")) {
            printf("Usage: ballsim_v2 [--config FILE] [--write-config [FILE]]\n");
            return 0;
        }
    }
    load_config(config_path);
    config_derive();

    srand((unsigned)time(NULL));

    cell_list_init(); // memset g_cell_head to -1 (done once here, not in loop)

#ifdef _OPENMP
    fprintf(stderr, "OpenMP: %d threads available.\n", omp_get_max_threads());
#endif

    if (SDL_Init(SDL_INIT_VIDEO) != 0) {
        fprintf(stderr, "SDL: %s\n", SDL_GetError());
        return 1;
    }
    SDL_Window   *win  = SDL_CreateWindow(
        "Parabolic Vessel \xe2\x80\x94 Hard Sphere Gas  [v2.0 \xe2\x80\x94 cell-list MD]",
        SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED,
        CFG_WIN_W, cfg.sim_h, SDL_WINDOW_SHOWN);
    SDL_Renderer *rend = SDL_CreateRenderer(
        win, -1, SDL_RENDERER_ACCELERATED | SDL_RENDERER_PRESENTVSYNC);
    if (!win || !rend) {
        fprintf(stderr, "SDL: %s\n", SDL_GetError());
        return 1;
    }

    int  n = cfg.n_particles;
    init_particles(n, true);
    sliders_init();

    Uint32 t_prev = SDL_GetTicks(), t_last_samp = t_prev;
    double accumulator = 0.0;
    const double FIXED_DT = 1.0 / 120.0;
    bool running = true;

    while (running) {
        SDL_Event ev;
        while (SDL_PollEvent(&ev)) {
            if (ev.type == SDL_QUIT) running = false;

            if (ev.type == SDL_MOUSEBUTTONDOWN && ev.button.button == SDL_BUTTON_LEFT) {
                int hit = slider_hit(ev.button.x, ev.button.y, g_slider_panel_top);
                if (hit >= 0) {
                    active_slider = hit;
                    if (slider_set(hit, ev.button.x, g_slider_panel_top)) {
                        n = (int)cfg.n_particles;
                        init_particles(n, true);
                    }
                }
            }
            if (ev.type == SDL_MOUSEBUTTONUP && ev.button.button == SDL_BUTTON_LEFT)
                active_slider = -1;
            if (ev.type == SDL_MOUSEMOTION && active_slider >= 0) {
                if (slider_set(active_slider, ev.motion.x, g_slider_panel_top)) {
                    n = (int)cfg.n_particles;
                    init_particles(n, true);
                }
            }

            if (ev.type == SDL_KEYDOWN)
                switch (ev.key.keysym.sym) {
                case SDLK_ESCAPE: case SDLK_q: running = false; break;
                case SDLK_r: init_particles(n, false); break;
                case SDLK_m: init_particles(n, true);  break;
                case SDLK_EQUALS: case SDLK_PLUS:
                    if (n < MAX_PARTICLES) { n++; cfg.n_particles = n; } break;
                case SDLK_MINUS:
                    if (n > 2) { n--; cfg.n_particles = n; } break;
                case SDLK_v:
                    if (!rec_active) { record_start(CFG_WIN_W,cfg.sim_h); rec_blink=SDL_GetTicks(); }
                    else record_stop();
                    break;
                case SDLK_c:
                    color_mode = (color_mode + 1) % N_COLOR_MODES;
                    break;
                }
        }

        Uint32 t_now = SDL_GetTicks();
        double frame_time = (t_now - t_prev) / 1000.0;
        if (frame_time > 0.05) frame_time = 0.05;
        t_prev = t_now;
        accumulator += frame_time;

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

        // Pre-compute normalisation for speed/KE color modes
        double max_val = 1.0;
        if (color_mode == 1) {
            for (int i = 0; i < n; i++) {
                double s = sqrt(p[i].vx*p[i].vx + p[i].vy*p[i].vy);
                if (s > max_val) max_val = s;
            }
        } else if (color_mode == 2) {
            for (int i = 0; i < n; i++) {
                double ke = 0.5*p[i].mass*(p[i].vx*p[i].vx + p[i].vy*p[i].vy);
                if (ke > max_val) max_val = ke;
            }
        }

        for (int i = 0; i < n; i++) {
            Uint8 cr, cg, cb;
            particle_rgb(i, max_val, &cr, &cg, &cb);
            draw_particle(rend, p[i].x, p[i].y, p[i].radius, cr, cg, cb);
        }
        draw_vessel(rend);
        double ke, pe;
        double te = compute_energy(n, &ke, &pe);
        draw_panel(rend, n, te, ke, pe);
        SDL_RenderPresent(rend);
        capture_frame(rend, CFG_WIN_W, cfg.sim_h);
    }

    record_stop();
    SDL_DestroyRenderer(rend);
    SDL_DestroyWindow(win);
    SDL_Quit();
    return 0;
}
