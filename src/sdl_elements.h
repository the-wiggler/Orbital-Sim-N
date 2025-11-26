#ifndef SDL_ELEMENTS_H
#define SDL_ELEMENTS_H

#include <math.h>
#include <SDL3/SDL.h>
#include <SDL3_ttf/SDL_ttf.h>
#include <stdbool.h>
#include "config.h"

extern TTF_Font* g_font;
extern TTF_Font* g_font_small;

typedef struct {
    int x, y, width, height;
    bool is_hovered;
    SDL_Color normal_color;
    SDL_Color hover_color;
} button_t;

typedef struct {
    button_t sc_button;
    button_t csv_load_button;
    button_t add_body_button;
    button_t add_craft_button;
    button_t show_stats_button;
} button_storage_t;

typedef enum {
    INPUT_NONE,
    INPUT_NAME,
    INPUT_MASS,
    INPUT_X_POS,
    INPUT_Y_POS,
    INPUT_X_VEL,
    INPUT_Y_VEL
} input_state_t;

typedef struct {
    bool active;
    input_state_t state;
    char input_buffer[256];
    char name[256];
    double mass;
    double x_pos;
    double y_pos;
    double x_vel;
    double y_vel;
} text_input_dialog_t;

void SDL_RenderFillCircle(SDL_Renderer* renderer, int centerX, int centerY, int radius);
void drawScaleBar(SDL_Renderer* renderer, window_params_t wp);
bool isMouseInRect(int mouse_x, int mouse_y, int rect_x, int rect_y, int rect_w, int rect_h);
void body_renderOrbitBodies(SDL_Renderer* renderer, body_properties_t* gb, int num_bodies, window_params_t wp);
void craft_renderCrafts(SDL_Renderer* renderer, spacecraft_properties_t* sc, int num_craft, window_params_t wp);
void renderTimeIndicators(SDL_Renderer* renderer, window_params_t wp);
void SDL_WriteText(SDL_Renderer* renderer, TTF_Font* font, const char* text, float x, float y, SDL_Color color);
void renderButton(SDL_Renderer* renderer, button_t* button, const char* text, window_params_t wp);
void renderUIButtons(SDL_Renderer* renderer, button_storage_t* buttons, window_params_t* wp);
void initButtons(button_storage_t* buttons, window_params_t wp);
void init_window_params(window_params_t* wp);
void init_text_dialog(text_input_dialog_t* dialog);
void displayError(const char* title, const char* message);
void showFPS(SDL_Renderer* renderer, Uint64 frame_start_timem, Uint64 perf_freq, window_params_t wp);

void runEventCheck(SDL_Event* event, window_params_t* wp, body_properties_t** bodies, int* num_bodies, spacecraft_properties_t** sc, int* num_craft, button_storage_t* buttons, text_input_dialog_t* dialog, stats_window_t* stats_window);
void renderBodyTextInputDialog(SDL_Renderer* renderer, text_input_dialog_t* dialog, window_params_t wp);
void renderStatsBox(SDL_Renderer* renderer, body_properties_t* bodies, int num_bodies, spacecraft_properties_t* sc, int num_craft, window_params_t wp, stats_window_t* stats_window);


#endif