#ifndef RENDERER_H
#define RENDERER_H

#include <SDL3/SDL.h>
#include "../types.h"

window_params_t init_window_params();
SDL_GL_init_t init_SDL_OPENGL_window(const char* title, int width, int height, Uint32* outWindowID);
void displayError(const char* title, const char* message);
void runEventCheck(SDL_Event* event, sim_properties_t* sim);
#endif
