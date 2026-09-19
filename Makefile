CC      = gcc
CFLAGS  = -O3 -std=c11 -Wall -Wextra -Wno-unused-parameter
SDL_CF  = $(shell pkg-config --cflags sdl2 2>/dev/null || echo -I/usr/include/SDL2)
SDL_LF  = $(shell pkg-config --libs   sdl2 2>/dev/null || echo -lSDL2)
LIBS    = $(SDL_LF) -lm

.PHONY: all clean run run2

all: ballsim ballsim_v2

# ── v1: event-driven, up to 500 particles ────────────────────────────────────
ballsim: verlet_collisions_two_gasses_video.c
	$(CC) $(CFLAGS) $(SDL_CF) -o $@ $< $(LIBS)
	@echo "Built v1: ./ballsim"

# ── v2: cell-list Verlet MD, up to 5000 particles, OpenMP ───────────────────
ballsim_v2: ballsim_v2.c
	$(CC) $(CFLAGS) $(SDL_CF) -fopenmp -o $@ $< $(LIBS)
	@echo "Built v2: ./ballsim_v2"

run:  ballsim   ; ./ballsim
run2: ballsim_v2; ./ballsim_v2

clean:
	rm -f ballsim ballsim_v2
