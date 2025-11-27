#include "config.h"
#ifdef _WIN32
    #include <windows.h>
#else
    #include <unistd.h>
#endif
#include <pthread.h>
#include <stdlib.h>
#include <SDL3/SDL.h>
#include <SDL3/SDL_main.h>
#include <SDL3_ttf/SDL_ttf.h>
#include <stdbool.h>
#include "sim_calculations.h"
#include "sdl_elements.h"

// NOTE: ALL CALCULATIONS SHOULD BE DONE IN BASE SI UNITS

// universal gravitation constant
char* FILENAME = "planet_data.json";
char* SPACECRAFT_FILENAME = "spacecraft_data.json";

TTF_Font* g_font = NULL;
TTF_Font* g_font_small = NULL;

pthread_mutex_t sim_vars_mutex;

////////////////////////////////////////////////////////////////////////////////////////////////////
// SIM CALCULATION FUNCTION
////////////////////////////////////////////////////////////////////////////////////////////////////
void* physicsSim(void* args) {
    const physics_sim_args* s = (physics_sim_args*)args;

    while (s->wp->sim_running) {
        // lock mutex before accessing data
        pthread_mutex_lock(&sim_vars_mutex);

        // IMPORTANT -- DOES ALL BODY CALCULATIONS:
        runCalculations(s->gb, s->sc, s->wp, *(s->num_bodies), *(s->num_craft));

        // unlock mutex when done :)
        pthread_mutex_unlock(&sim_vars_mutex);
    }

    return NULL;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// MAIN
////////////////////////////////////////////////////////////////////////////////////////////////////
int main() {
    ////////////////////////////////////////
    // INIT                               //
    ////////////////////////////////////////
    // initialize simulation parameters
    window_params_t wp = {0};
    init_window_params(&wp);

    // initialize UI elements
    button_storage_t buttons;
    initButtons(&buttons, wp);

    // initialize text input
    text_input_dialog_t dialog = {0};
    init_text_dialog(&dialog);

    // initialize SDL3
    SDL_Init(SDL_INIT_VIDEO);
    // create an SDL window
    SDL_Window* window = SDL_CreateWindow("Orbit Simulation", (int)wp.window_size_x, (int)wp.window_size_y, SDL_WINDOW_RESIZABLE);
    wp.main_window_ID = SDL_GetWindowID(window);
    // create an SDL renderer and clear the window to create a blank canvas
    SDL_Renderer *renderer = SDL_CreateRenderer(window, NULL);
    SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
    SDL_RenderClear(renderer);

    // fps counter init
    Uint64 perf_freq = SDL_GetPerformanceFrequency();
    Uint64 frame_start;

    // SDL ttf font stuff
    TTF_Init();
    g_font = TTF_OpenFont("CascadiaCode.ttf", wp.font_size);
    g_font_small = TTF_OpenFont("CascadiaCode.ttf", (float)wp.window_size_x / 90);

    // toggleable stats window
    stats_window_t stats_window = {0};

    ////////////////////////////////////////
    // SIM VARS                           //
    ////////////////////////////////////////

    // initialize simulation objects
    int num_bodies = 0;
    body_properties_t* gb = NULL;
    int num_craft = 0;
    spacecraft_properties_t* sc = NULL;

    ////////////////////////////////////////
    // SIM THREAD INIT                    //
    ////////////////////////////////////////
    // initialize simulation thread
    pthread_t simThread;
    pthread_mutex_init(&sim_vars_mutex, NULL);
    
    // arguments to pass into the sim thread
    physics_sim_args ps_args = {
        .gb = &gb,
        .sc = &sc,
        .wp = &wp,
        .num_bodies = &num_bodies,
        .num_craft = &num_craft
    };

    // creates the sim thread
    if (pthread_create(&simThread, NULL, physicsSim, &ps_args) != 0) {
        displayError("ERROR", "Error when creating physics simulation process");
        return 1;
    }

    ////////////////////////////////////////////////////////
    // simulation loop                                    //
    ////////////////////////////////////////////////////////
    while (wp.window_open) {
        // measure frame start time
        frame_start = SDL_GetPerformanceCounter();

        // clears previous frame from the screen
        SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
        SDL_RenderClear(renderer);

        // draw scale reference bar
        drawScaleBar(renderer, wp);

        // draw speed control button
        renderUIButtons(renderer, &buttons, &wp);

        // render text input dialog if active
        renderBodyTextInputDialog(renderer, &dialog, wp);

        // draw time indicator text
        renderTimeIndicators(renderer, wp);

        // lock body_sim mutex for reading
        pthread_mutex_lock(&sim_vars_mutex);

        // user input event checking logic
        SDL_Event event;
        runEventCheck(&event, &wp, &gb, &num_bodies, &sc, &num_craft, &buttons, &dialog, &stats_window);

        // render the bodies
        body_renderOrbitBodies(renderer, gb, num_bodies, wp);

        // render the spacecraft
        craft_renderCrafts(renderer, sc, num_craft);

        // render stats in main window if enabled
        if (stats_window.is_shown) {
            renderStatsBox(renderer, gb, num_bodies, sc, num_craft, wp, &stats_window);
        }

        // unlock sim vars mutex when done
        pthread_mutex_unlock(&sim_vars_mutex);

        // shows and limits FPS
        showFPS(renderer, frame_start, perf_freq, wp);

        // present the renderer to the screen
        SDL_RenderPresent(renderer);
    }
    ////////////////////////////////////////////////////////
    // end of simulation loop                             //
    ////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////
    // CLEAN UP                                       //
    ////////////////////////////////////////////////////
    // wait for simulation thread to finish
    pthread_join(simThread, NULL);

    // destroy mutex
    pthread_mutex_destroy(&sim_vars_mutex);

    if (gb != NULL) {
        for (int i = 0; i < num_bodies; i++) {
            free(gb[i].name);
        }
        free(gb);
        gb = NULL;
    }
    num_bodies = 0;

    if (sc != NULL) {
        for (int i = 0; i < num_craft; i++) {
            free(sc[i].name);
        }
        free(sc);
        sc = NULL;
    }
    num_craft = 0;

    if (g_font) TTF_CloseFont(g_font);
    if (g_font_small) TTF_CloseFont(g_font_small);
    TTF_Quit();
    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();
    return 0;
}