#ifndef GLOBALS_H
#define GLOBALS_H

#include <SDL3/SDL.h>
static const double G = 6.67430E-11;

static const double PI = 3.14159265358979323846;
static const float PI_f = 3.14159265358979323846f;

static const double PI_OVER_6    =  0.52359877559829877312;  // π/6
static const double PI_OVER_4    =  0.78539816339744830962;  // π/4
static const double PI_OVER_3    =  1.04719755119659763133;  // π/3
static const double PI_OVER_2    =  1.57079632679489661923;  // π/2
static const double TWO_PI       =  6.28318530717958623200;  // 2π

static const float PI_OVER_6_f   =  0.5235987756f;  // π/6
static const float PI_OVER_4_f   =  0.7853981634f;  // π/4
static const float PI_OVER_3_f   =  1.0471975512f;  // π/3
static const float PI_OVER_2_f   =  1.57079632679f; // π/2
static const float TWO_PI_f      =  6.28318530718f; // 2π

static const char* SIMULATION_FILENAME = "simulation_data.json";
static const float SCALE = 1e7f; // scales in-sim meters to openGL coordinates -- this is an arbitrary number that can be adjusted
static const int MAX_PLANETS = 16;
static const int PATH_CAPACITY = 1000;

static const SDL_Color TEXT_COLOR = {210, 210, 210, 255};
static const SDL_Color BUTTON_COLOR = {30,30,30, 255};
static const SDL_Color BUTTON_HOVER_COLOR = {20,20,20, 255};
static const SDL_Color ACCENT_COLOR = {80, 150, 220, 255};


#endif
