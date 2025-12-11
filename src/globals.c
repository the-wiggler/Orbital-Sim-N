#include "globals.h"

TTF_Font* g_font = NULL;
TTF_Font* g_font_small = NULL;
pthread_mutex_t sim_vars_mutex;
const double G = 6.67430E-11;
char* SIMULATION_FILENAME = "simulation_data.json";

SDL_Color TEXT_COLOR = {210, 210, 210, 255};
SDL_Color BUTTON_COLOR = {70,70,70, 255};
SDL_Color BUTTON_HOVER_COLOR = {50,50,50, 255};
SDL_Color ACCENT_COLOR = {80, 150, 220, 255};
