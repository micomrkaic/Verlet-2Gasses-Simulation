CC     = gcc
TARGET = ballsim

UNAME := $(shell uname -s)
ifeq ($(UNAME), Darwin)
    SDL_CFLAGS = -I/opt/homebrew/include/SDL2 -I/usr/local/include/SDL2
    SDL_LIBS   = -L/opt/homebrew/lib -L/usr/local/lib -lSDL2
else
    SDL_CFLAGS = $(shell pkg-config --cflags sdl2 2>/dev/null || echo -I/usr/include/SDL2)
    SDL_LIBS   = $(shell pkg-config --libs   sdl2 2>/dev/null || echo -lSDL2)
endif

CFLAGS = -O2 -std=c11 -Wall -Wextra -Wno-unused-parameter $(SDL_CFLAGS)
LIBS   = $(SDL_LIBS) -lm

.PHONY: all run clean

all: $(TARGET)

$(TARGET): verlet_collisions_two_gasses_video.c
	$(CC) $(CFLAGS) -o $@ $< $(LIBS)
	@echo "Built: ./$(TARGET)"

run: all
	./$(TARGET)

clean:
	rm -f $(TARGET)
