#ifndef SDL_ELEMENTS_H
#define SDL_ELEMENTS_H

#include <SDL3/SDL.h>
#include <SDL3_ttf/SDL_ttf.h>
#include <stdbool.h>
#include "config.h"

extern TTF_Font* g_font;
extern TTF_Font* g_font_small;

typedef struct {
    float x, y, width, height;
    bool is_hovered;
    SDL_Color normal_color;
    SDL_Color hover_color;
} button_t;

typedef struct {
    button_t sc_button;
    button_t csv_load_button;
    button_t add_craft_button;
    button_t show_stats_button;
} button_storage_t;

void SDL_RenderFillCircle(SDL_Renderer* renderer, float centerX, float centerY, float radius);
void drawScaleBar(SDL_Renderer* renderer, window_params_t wp);
bool isMouseInRect(int mouse_x, int mouse_y, int rect_x, int rect_y, int rect_w, int rect_h);
void body_renderOrbitBodies(SDL_Renderer* renderer, body_properties_t* gb, window_params_t wp);
void craft_renderCrafts(SDL_Renderer* renderer, const spacecraft_properties_t* sc);
void renderTimeIndicators(SDL_Renderer* renderer, window_params_t wp);
void SDL_WriteText(SDL_Renderer* renderer, TTF_Font* font, const char* text, float x, float y, SDL_Color color);
void renderButton(SDL_Renderer* renderer, const button_t* button, const char* text);
void renderUIButtons(SDL_Renderer* renderer, const button_storage_t* buttons, const window_params_t* wp);
void initButtons(button_storage_t* buttons, window_params_t wp);
void init_window_params(window_params_t* wp);
void displayError(const char* title, const char* message);
void showFPS(SDL_Renderer* renderer, Uint64 frame_start_time, Uint64 perf_freq, window_params_t wp, bool FPS_SHOWN);

void runEventCheck(SDL_Event* event, window_params_t* wp, body_properties_t* bodies, spacecraft_properties_t* sc, button_storage_t* buttons, stats_window_t* stats_window);
void renderStatsBox(SDL_Renderer* renderer, body_properties_t* bodies, const spacecraft_properties_t* sc, window_params_t wp, stats_window_t* stats_window);


#endif