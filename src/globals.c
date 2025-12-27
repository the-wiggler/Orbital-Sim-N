#include "globals.h"

pthread_mutex_t sim_vars_mutex;
const double G = 6.67430E-11;
const double M_PI = 3.14159265358979323846;
char* SIMULATION_FILENAME = "simulation_data.json";

float SCALE = 1e7; // scales in-sim meters to openGL coordinates -- this is an arbitrary number that can be adjusted

SDL_Color TEXT_COLOR = {210, 210, 210, 255};
SDL_Color BUTTON_COLOR = {30,30,30, 255};
SDL_Color BUTTON_HOVER_COLOR = {20,20,20, 255};
SDL_Color ACCENT_COLOR = {80, 150, 220, 255};
