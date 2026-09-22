CC      = gcc
CFLAGS  = -O3 -std=c11 -Wall -Wextra -Wno-unused-parameter
UNAME  := $(shell uname -s)

ifeq ($(UNAME), Darwin)
    SDL_CF = -I/opt/homebrew/include -I/opt/homebrew/include/SDL2 \
             -I/usr/local/include   -I/usr/local/include/SDL2
    SDL_LF = -L/opt/homebrew/lib -L/usr/local/lib -lSDL2

    # Match only gcc-NN (digits only after the dash) — avoids gcc-ranlib-NN etc.
    GCC_BIN := $(shell find /opt/homebrew/bin /usr/local/bin -name 'gcc-[0-9]*' \
                    -not -name 'gcc-*-*' 2>/dev/null | sort | tail -1)

    ifneq ($(GCC_BIN),)
        CC_V2  = $(GCC_BIN)
        OMP_CF = -fopenmp
        OMP_LF = -fopenmp
    else
        CC_V2   = $(CC)
        LIBOMP := $(shell brew --prefix libomp 2>/dev/null)
        ifneq ($(LIBOMP),)
            OMP_CF = -Xpreprocessor -fopenmp -I$(LIBOMP)/include
            OMP_LF = -L$(LIBOMP)/lib -lomp
        else
            $(warning OpenMP not found — ballsim_v2 will be single-threaded.)
            $(warning Fix: brew install gcc   OR   brew install libomp)
            OMP_CF =
            OMP_LF =
        endif
    endif
else
    SDL_CF = $(shell pkg-config --cflags sdl2 2>/dev/null || echo -I/usr/include/SDL2)
    SDL_LF = $(shell pkg-config --libs   sdl2 2>/dev/null || echo -lSDL2)
    CC_V2  = $(CC)
    OMP_CF = -fopenmp
    OMP_LF = -fopenmp
endif

LIBS = $(SDL_LF) -lm

.PHONY: all clean run run2

all: ballsim ballsim_v2

ballsim: verlet_collisions_two_gasses_video.c
	$(CC) $(CFLAGS) $(SDL_CF) -o $@ $< $(LIBS)
	@echo "Built v1: ./ballsim"

ballsim_v2: ballsim_v2.c
	$(CC_V2) $(CFLAGS) $(SDL_CF) $(OMP_CF) -o $@ $< $(OMP_LF) $(LIBS)
	@echo "Built v2: ./ballsim_v2"

run:  ballsim   ; ./ballsim
run2: ballsim_v2; ./ballsim_v2

clean:
	rm -f ballsim ballsim_v2
