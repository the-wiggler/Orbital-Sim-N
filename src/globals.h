#ifndef GLOBALS_H
#define GLOBALS_H

#include <SDL3/SDL.h>
#define G 6.67430E-11
#define PI 3.14159265358979323846
#define M_PI_f 3.14159265358979323846f
static const char* SIMULATION_FILENAME = "simulation_data.json";
#define SCALE 1e7f // scales in-sim meters to openGL coordinates -- this is an arbitrary number that can be adjusted

static const SDL_Color TEXT_COLOR = {210, 210, 210, 255};
static const SDL_Color BUTTON_COLOR = {30,30,30, 255};
static const SDL_Color BUTTON_HOVER_COLOR = {20,20,20, 255};
static const SDL_Color ACCENT_COLOR = {80, 150, 220, 255};


#endif
