#ifndef RENDERER_H
#define RENDERER_H

#include <SDL3/SDL.h>
#include "../types.h"

window_params_t init_window_params(void);
console_t init_console(window_params_t window_params);
SDL_GL_init_t init_SDL_OPENGL_window(const char* title, sim_properties_t* sim);
void displayError(const char* title, const char* message);
void runEventCheck(SDL_Event* event, sim_properties_t* sim, binary_filenames_t* filenames);
void renderCMDWindow(sim_properties_t sim, font_t* font);

// GL scene colors
typedef struct { float r, g, b; } gl_color3_t;
typedef struct { float r, g, b, a; } gl_color4_t;

// background
static const gl_color3_t BG_COLOR = {0.0F, 0.0F, 0.0F};

// coordinate grid
static const gl_color3_t GRID_X_AXIS_COLOR = {0.4F, 0.4F, 0.4F};
static const gl_color3_t GRID_Y_AXIS_COLOR = {0.4F, 0.4F, 0.4F};
static const gl_color3_t GRID_Z_AXIS_COLOR = {0.4F, 0.4F, 0.4F};
static const gl_color3_t GRID_DIAGONAL_COLOR = {0.2F, 0.2F, 0.2F};

// sphere of influence overlay
static const gl_color4_t SOI_OVERLAY_COLOR = {0.0F, 0.5F, 1.0F, 0.2F};

// orbit trails
static const gl_color3_t BODY_ORBIT_COLOR = {0.0F, 0.6F, 0.6F};
static const gl_color3_t CRAFT_ORBIT_COLOR = {0.0F, 1.0F, 1.0F};
static const gl_color3_t CRAFT_PATH_COLOR = {1.0F, 1.0F, 1.0F};

// debug visuals
static const gl_color3_t BODY_LINE_COLOR = {0.0F, 1.0F, 1.0F};
static const gl_color3_t INCLINATION_ABOVE_COLOR = {0.5F, 0.5F, 1.0F};
static const gl_color3_t INCLINATION_BELOW_COLOR = {1.0F, 0.5F, 0.5F};
static const gl_color3_t ROTATION_AXIS_COLOR = {0.0F, 1.0F, 1.0F};
static const gl_color3_t CRAFT_RADIUS_COLOR = {1.0F, 1.0F, 1.0F};
static const gl_color3_t CRAFT_EQUATORIAL_COLOR = {1.0F, 0.0F, 0.0F};
static const gl_color3_t CRAFT_INCLINATION_COLOR = {0.0F, 1.0F, 0.0F};

// HUD text
static const gl_color3_t HUD_TEXT_COLOR = {1.0F, 1.0F, 1.0F};

#endif
