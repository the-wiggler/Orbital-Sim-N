#ifndef RENDERER_H
#define RENDERER_H

#include <SDL3/SDL.h>
#include "../types.h"

window_params_t init_window_params(void);
console_t init_console(window_params_t wp);
SDL_GL_init_t init_SDL_OPENGL_window(const char* title, int width, int height, Uint32* outWindowID);
void displayError(const char* title, const char* message);
void runEventCheck(SDL_Event* event, sim_properties_t* sim);
void renderCMDWindow(sim_properties_t* sim, font_t* font);

static const SDL_Color TEXT_COLOR = {210, 210, 210, 255};
static const SDL_Color BUTTON_COLOR = {30,30,30, 255};
static const SDL_Color BUTTON_HOVER_COLOR = {20,20,20, 255};
static const SDL_Color ACCENT_COLOR = {80, 150, 220, 255};

#endif
