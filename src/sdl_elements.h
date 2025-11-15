#ifndef SDL_ELEMENTS_H
#define SDL_ELEMENTS_H

#include <stdio.h>
#include <math.h>
#include <SDL3/SDL.h>
#include <SDL3_ttf/SDL_ttf.h>
#include <string.h>
#include <stdbool.h>
#include <stdlib.h>
#include "config.h"

extern TTF_Font* g_font;

typedef struct {
    int x, y, width, height;
    bool is_hovered;
} speed_control_t;

void SDL_RenderFillCircle(SDL_Renderer* renderer, int centerX, int centerY, int radius);
void drawScaleBar(SDL_Renderer* renderer, window_params_t wp);
bool isMouseInRect(int mouse_x, int mouse_y, int rect_x, int rect_y, int rect_w, int rect_h);
void drawSpeedControl(SDL_Renderer* renderer, speed_control_t* control, window_params_t wp);
void drawStatsBox(SDL_Renderer* renderer, body_properties_t* bodies, int num_bodies, double sim_time, window_params_t wp);
void SDL_WriteText(SDL_Renderer* renderer, TTF_Font* font, const char* text, float x, float y, SDL_Color color);

void runEventCheck(SDL_Event* event, speed_control_t* speed_control, window_params_t* wp, body_properties_t** bodies, int* num_bodies);


#endif