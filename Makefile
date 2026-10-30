CC      = gcc
CFLAGS  = -O3 -std=c11 -Wall -Wextra -Wno-unused-parameter
UNAME  := $(shell uname -s)
VERSION = 2.2

# ── Install locations (XDG user layout; override: make install PREFIX=/usr/local)
PREFIX  ?= $(HOME)/.local
BINDIR  ?= $(PREFIX)/bin
DATADIR ?= $(PREFIX)/share/ballsim
MANDIR  ?= $(PREFIX)/share/man/man1
INSTALL ?= install

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
# Compile the install-time asset directory into v2 so an installed binary
# finds its default config from any working directory.
DEFS = -DDATADIR='"$(DATADIR)"' -DBALLSIM_VERSION='"$(VERSION)"'

.PHONY: all clean run run2 bench install uninstall

all: ballsim ballsim_v2

# ── v1: event-driven, exact CCD, up to 500 particles ─────────────────────────
ballsim: verlet_collisions_two_gasses_video.c
	$(CC) $(CFLAGS) $(SDL_CF) -o $@ $< $(LIBS)
	@echo "Built v1: ./ballsim"

# ── v2: cell-list Verlet MD, up to 5000 particles, OpenMP ───────────────────
ballsim_v2: ballsim_v2.c
	$(CC_V2) $(CFLAGS) $(SDL_CF) $(OMP_CF) $(DEFS) -o $@ $< $(OMP_LF) $(LIBS)
	@echo "Built v2: ./ballsim_v2"

run:  ballsim    ; ./ballsim
run2: ballsim_v2 ; ./ballsim_v2

# Headless benchmark sweep (no window). Override: make bench STEPS=2400
STEPS ?= 1200
bench: ballsim_v2
	@for n in 300 1000 2000 5000; do ./ballsim_v2 --bench $$n $(STEPS) --seed 1; done

# ── Install ──────────────────────────────────────────────────────────────────
#   $(BINDIR)/ballsim, ballsim_v2          binaries
#   $(DATADIR)/ballsim_v2.cfg              pristine default config (generated)
#   $(MANDIR)/ballsim_v2.1, ballsim.1      man page (+ alias)
# Personal edits belong in ~/.config/ballsim/ballsim_v2.cfg, which the binary
# prefers over the share copy — see ballsim_v2(1), FILES.
install: all
	$(INSTALL) -d "$(BINDIR)" "$(DATADIR)" "$(MANDIR)"
	$(INSTALL) -m 755 ballsim ballsim_v2 "$(BINDIR)/"
	./ballsim_v2 --write-config "$(DATADIR)/ballsim_v2.cfg"
	$(INSTALL) -m 644 ballsim_v2.1 "$(MANDIR)/ballsim_v2.1"
	printf '.so man1/ballsim_v2.1\n' > "$(MANDIR)/ballsim.1"
	@echo
	@echo "Installed to $(PREFIX):"
	@echo "  $(BINDIR)/ballsim_v2   $(BINDIR)/ballsim"
	@echo "  $(DATADIR)/ballsim_v2.cfg"
	@echo "  $(MANDIR)/ballsim_v2.1"
	@case ":$$PATH:" in *":$(BINDIR):"*) ;; *) \
	  echo; echo "NOTE: $(BINDIR) is not in PATH. Add to your shell rc:"; \
	  echo "  export PATH=\"$(BINDIR):\$$PATH\"";; esac
	@echo "If 'man ballsim_v2' is not found: export MANPATH=\"$(PREFIX)/share/man:\$$MANPATH\""

uninstall:
	rm -f "$(BINDIR)/ballsim" "$(BINDIR)/ballsim_v2"
	rm -f "$(DATADIR)/ballsim_v2.cfg"
	rm -f "$(MANDIR)/ballsim_v2.1" "$(MANDIR)/ballsim.1"
	-rmdir "$(DATADIR)" 2>/dev/null
	@echo "Uninstalled from $(PREFIX). Your ~/.config/ballsim/ was left untouched."

clean:
	rm -f ballsim ballsim_v2
