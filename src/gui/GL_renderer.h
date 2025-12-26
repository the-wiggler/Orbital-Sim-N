//
// Created by java on 12/25/25.
//

#ifndef ORBITSIMULATION_GL_RENDERER_H
#define ORBITSIMULATION_GL_RENDERER_H

#include <SDL3/SDL.h>
#include <GL/glew.h>
#include <GL/gl.h>

static void update_viewport(SDL_Window *window);
char* loadShaderSource(const char* filepath);

// Transform world coordinates to normalized device coordinates (NDC)
void worldToNDC(double world_x, double world_y,
                double camera_x, double camera_y,
                double meters_per_pixel,
                float screen_width, float screen_height,
                float* out_x, float* out_y);

#endif //ORBITSIMULATION_GL_RENDERER_H