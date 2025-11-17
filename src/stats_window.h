#ifndef STATS_WINDOW_H
#define STATS_WINDOW_H

#include "config.h"
#include <SDL3/SDL.h>
#include <SDL3_ttf/SDL_ttf.h>

bool statsWindowInit(stats_window_t* stats);
void StatsWindow_handleEvent(stats_window_t* stats, SDL_Event* e);
void StatsWindow_render(stats_window_t* stats, int fps, float posX, float posY, body_properties_t* gb, int num_bodies, window_params_t wp);
void StatsWindow_show(stats_window_t* stats);
void StatsWindow_destroy(stats_window_t* stats);


#endif