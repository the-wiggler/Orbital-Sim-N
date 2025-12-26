#ifndef RENDERER_H
#define RENDERER_H

#include <SDL3/SDL.h>
#include "../types.h"

void init_window_params(window_params_t* wp);
void displayError(const char* title, const char* message);
void runEventCheck(SDL_Event* event, sim_properties_t* sim);
#endif
