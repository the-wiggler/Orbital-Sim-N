#ifndef GLOBALS_H
#define GLOBALS_H

#include <pthread.h>
#include <SDL3/SDL.h>

extern pthread_mutex_t sim_vars_mutex;
extern const double G;
extern const double M_PI;
extern char* SIMULATION_FILENAME;

extern SDL_Color ACCENT_COLOR;

#endif
