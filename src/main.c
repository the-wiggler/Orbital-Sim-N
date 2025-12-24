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
    sim_properties_t* sim = (sim_properties_t*)args;
    while (sim->wp.window_open) {
        while (sim->wp.sim_running) {
            // lock mutex before accessing data
            pthread_mutex_lock(&sim_vars_mutex);

            // DOES ALL BODY AND CRAFT CALCULATIONS:
            runCalculations(sim);

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
    ///
    // initialize simulation objects
    sim_properties_t sim = {
        .gb = {0},
        .gs = {0},
        .wp = {0}
    };

    // binary file creation
    binary_filenames_t filenames = {
        .global_data_FILE = fopen("global_data.bin", "wb")
    };

    // initialize window parameters
    SDL_Init(SDL_INIT_VIDEO);
    init_window_params(&sim.wp);

    // initialize UI elements
    button_storage_t buttons;
    initButtons(&buttons, sim.wp);
    stats_window_t stats_window = {0};

    // initialize SDL3 window
    // create an SDL window
    SDL_Window* window = SDL_CreateWindow("Orbit Simulation", (int)sim.wp.window_size_x, (int)sim.wp.window_size_y, SDL_WINDOW_RESIZABLE);
    sim.wp.main_window_ID = SDL_GetWindowID(window);
    // create an SDL renderer and clear the window to create a blank canvas
    SDL_Renderer *renderer = SDL_CreateRenderer(window, NULL);
    SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
    SDL_RenderClear(renderer);

    // fps counter init
    Uint64 perf_freq = SDL_GetPerformanceFrequency();
    Uint64 frame_start;

    // SDL ttf font stuff
    TTF_Init();
    g_font = TTF_OpenFont("font.ttf", sim.wp.font_size);
    g_font_small = TTF_OpenFont("font.ttf", (float)sim.wp.window_size_x / 90);

    ////////////////////////////////////////
    // SIM THREAD INIT                    //
    ////////////////////////////////////////
    // initialize simulation thread
    pthread_t simThread;
    pthread_mutex_init(&sim_vars_mutex, NULL);

    // creates the sim thread
    if (pthread_create(&simThread, NULL, physicsSim, &sim) != 0) {
        displayError("ERROR", "Error when creating physics simulation process");
        sim.wp.sim_running = false;
        sim.wp.window_open = false;
        return 1;
    }

    ////////////////////////////////////////////////////////
    // simulation loop                                    //
    ////////////////////////////////////////////////////////
    while (sim.wp.window_open) {

        // measure frame start time
        frame_start = SDL_GetPerformanceCounter();

        // clears previous frame from the screen
        SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
        SDL_RenderClear(renderer);

        // lock body_sim mutex for reading
        pthread_mutex_lock(&sim_vars_mutex);

        // user input event checking logic
        SDL_Event event;
        runEventCheck(&event, &sim, &buttons, &stats_window);

        // draw speed control button
        renderUIButtons(renderer, &buttons, &sim.wp);

        // draw time indicator text
        renderTimeIndicators(renderer, sim.wp);

        // if the main view is shown, do the respective rendering functions
        if (sim.wp.main_view_shown && !sim.wp.craft_view_shown) {
            // draw scale reference bar
            drawScaleBar(renderer, sim.wp);

            // render the bodies
            body_renderOrbitBodies(renderer, &sim);
//
            // render the spacecraft
            craft_renderCrafts(renderer, &sim.gs);

            // render stats in main window if enabled
            if (stats_window.is_shown) renderStatsBox(renderer, &sim, &stats_window);
        }

        // if the craft view is shown, do the respective rendering functions
        if (sim.wp.craft_view_shown) {
            craft_RenderCraftView(renderer, &sim);
        }

        if (sim.wp.data_logging_enabled) {
            exportTelemetryBinary(filenames, &sim);
        }

        // check if sim needs to be reset
        if (sim.wp.reset_sim) resetSim(&sim);

        // unlock sim vars mutex when done
        pthread_mutex_unlock(&sim_vars_mutex);

        // limits FPS
        showFPS(renderer, frame_start, perf_freq, sim.wp, false);

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
    cleanup(&sim);

    // cleanup button textures
    destroyAllButtonTextures(&buttons);

    fclose(filenames.global_data_FILE);
    if (g_font) TTF_CloseFont(g_font);
    if (g_font_small) TTF_CloseFont(g_font_small);
    TTF_Quit();
    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();
    return 0;
}