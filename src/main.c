#include <stdio.h>

#include "globals.h"
#include "types.h"
#include "sim/simulation.h"
#include "gui/renderer.h"
#include "gui/craft_view.h"
#include "utility/telemetry_export.h"
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

// NOTE: ALL CALCULATIONS SHOULD BE DONE IN BASE SI UNITS

////////////////////////////////////////////////////////////////////////////////////////////////////
// PHYSICS SIMULATION THREAD
////////////////////////////////////////////////////////////////////////////////////////////////////
void* physicsSim(void* args) {
    const physics_sim_args* s = (physics_sim_args*)args;
    while (s->wp->window_open) {
        while (s->wp->sim_running) {
            // lock mutex before accessing data
            pthread_mutex_lock(&sim_vars_mutex);

            // DOES ALL BODY AND CRAFT CALCULATIONS:
            runCalculations(s->gb, s->sc, s->wp);

            // unlock mutex when done :)
            pthread_mutex_unlock(&sim_vars_mutex);
        }
    }
    return NULL;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// MAIN :)
////////////////////////////////////////////////////////////////////////////////////////////////////
int main(int argc, char *argv[]) {
    (void)argc;
    (void)argv;
    ////////////////////////////////////////
    // INIT                               //
    ////////////////////////////////////////

    // binary file creation
    binary_filenames_t filenames = {
        .body_pos_FILE = fopen("body_pos_data.bin", "wb"),
    };
    if (!filenames.body_pos_FILE) displayError("ERROR", "Error: unable to create position data binary file");

    // initialize window parameters
    SDL_Init(SDL_INIT_VIDEO);
    window_params_t wp = {0};
    init_window_params(&wp);

    // initialize UI elements
    button_storage_t buttons;
    initButtons(&buttons, wp);
    stats_window_t stats_window = {0};

    // initialize SDL3 window
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
    g_font = TTF_OpenFont("font.ttf", wp.font_size);
    g_font_small = TTF_OpenFont("font.ttf", (float)wp.window_size_x / 90);

    ////////////////////////////////////////
    // SIM VARS                           //
    ////////////////////////////////////////

    // initialize simulation objects
    body_properties_t gb = {0};
    spacecraft_properties_t sc = {0};

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
        .wp = &wp
    };

    // creates the sim thread
    if (pthread_create(&simThread, NULL, physicsSim, &ps_args) != 0) {
        displayError("ERROR", "Error when creating physics simulation process");
        wp.sim_running = false;
        wp.window_open = false;
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

        // lock body_sim mutex for reading
        pthread_mutex_lock(&sim_vars_mutex);

        // user input event checking logic
        SDL_Event event;
        runEventCheck(&event, &wp, &gb, &sc, &buttons, &stats_window);

        // draw speed control button
        renderUIButtons(renderer, &buttons, &wp);

        // draw time indicator text
        renderTimeIndicators(renderer, wp);

        // if the main view is shown, do the respective rendering functions
        if (wp.main_view_shown && !wp.craft_view_shown) {
            // draw scale reference bar
            drawScaleBar(renderer, wp);

            // render the bodies
            body_renderOrbitBodies(renderer, &gb, &wp);
//
            // render the spacecraft
            craft_renderCrafts(renderer, &sc);

            // render stats in main window if enabled
            if (stats_window.is_shown) renderStatsBox(renderer, &gb, &sc, wp, &stats_window);
        }

        // if the craft view is shown, do the respective rendering functions
        if (wp.craft_view_shown) {
            craft_RenderCraftView(renderer, &wp, &gb, &sc);
        }

        wp.data_logging_enabled = true;
        if (wp.data_logging_enabled) {
            exportTelemetryBinary(filenames, wp, gb, sc);
        }

        // check if sim needs to be reset
        if (wp.reset_sim) resetSim(&wp, &gb, &sc);

        // unlock sim vars mutex when done
        pthread_mutex_unlock(&sim_vars_mutex);

        // limits FPS
        showFPS(renderer, frame_start, perf_freq, wp, false);

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

    // cleanup all allocated memory
    cleanup_args clean_args = {
        .gb = &gb,
        .sc = &sc,
        .wp = &wp,
    };
    cleanup(&clean_args);

    // cleanup button textures
    destroyAllButtonTextures(&buttons);

    fclose(filenames.body_pos_FILE);
    if (g_font) TTF_CloseFont(g_font);
    if (g_font_small) TTF_CloseFont(g_font_small);
    TTF_Quit();
    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();
    return 0;
}