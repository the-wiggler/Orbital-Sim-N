#include "stats_window.h"
#include "config.h"
#include "sdl_elements.h"

// initialize the stats window
bool statsWindowInit(stats_window_t* stats) {
    if (!stats) return false;

    stats->window = SDL_CreateWindow("Stats Window", 400, 300, 0);
    if (!stats->window) {
        SDL_Log("Failed to create stats window: %s", SDL_GetError());
        return false;
    }

    stats->renderer = SDL_CreateRenderer(stats->window, NULL);
    if (!stats->renderer) {
        SDL_Log("Failed to create renderer: %s", SDL_GetError());
        SDL_DestroyWindow(stats->window);
        return false;
    }

    stats->window_ID = SDL_GetWindowID(stats->window);
    stats->is_shown = true;
    return true;
}

// handle events for the stats window
void StatsWindow_handleEvent(stats_window_t* stats, SDL_Event* e) {
    // check if event belongs to this window
    if (e->type >= SDL_EVENT_WINDOW_FIRST && 
        e->type <= SDL_EVENT_WINDOW_LAST &&
        e->window.windowID == stats->window_ID) {
        
        if (e->type == SDL_EVENT_WINDOW_CLOSE_REQUESTED) {
            SDL_HideWindow(stats->window);
            stats->is_shown = false;
        }
    }
}

// render the stats window with data from main window
void StatsWindow_render(stats_window_t* stats, int fps, float posX, float posY, body_properties_t* gb, int num_bodies, window_params_t wp) {
    if (!stats->is_shown) return;
    
    SDL_SetRenderDrawColor(stats->renderer, 30, 30, 30, 255);
    SDL_RenderClear(stats->renderer);

    // draw the stats
    drawStatsBox(stats->renderer, gb, num_bodies, wp.sim_time, wp);
    
    SDL_RenderPresent(stats->renderer);
}

// show the window
void StatsWindow_show(stats_window_t* stats) {
    SDL_ShowWindow(stats->window);
    stats->is_shown = true;
}

// cleanup
void StatsWindow_destroy(stats_window_t* stats) {
    if (stats->renderer) {
        SDL_DestroyRenderer(stats->renderer);
    }
    if (stats->window) {
        SDL_DestroyWindow(stats->window);
    }
}