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
    while (s->wp->window_open) {
        while (s->wp->sim_running) {
            // lock mutex before accessing data
            pthread_mutex_lock(&sim_vars_mutex);

            // IMPORTANT -- DOES ALL BODY CALCULATIONS:
            runCalculations(s->gb, s->sc, s->wp);

            // unlock mutex when done :)
            pthread_mutex_unlock(&sim_vars_mutex);
        }
    }

    return NULL;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// MAIN
////////////////////////////////////////////////////////////////////////////////////////////////////
int main(int argc, char *argv[]) {
    (void)argc;  // Unused parameter
    (void)argv;  // Unused parameter
    ////////////////////////////////////////
    // INIT                               //
    ////////////////////////////////////////

    // initialize window parameters
    SDL_Init(SDL_INIT_VIDEO);
    window_params_t wp = {0};
    init_window_params(&wp);

    // initialize UI elements
    button_storage_t buttons;
    initButtons(&buttons, wp);

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
    g_font = TTF_OpenFont("CascadiaCode.ttf", wp.font_size);
    g_font_small = TTF_OpenFont("CascadiaCode.ttf", (float)wp.window_size_x / 90);

    // toggleable stats window
    stats_window_t stats_window = {0};

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

        // draw time indicator text
        renderTimeIndicators(renderer, wp);

        // lock body_sim mutex for reading
        pthread_mutex_lock(&sim_vars_mutex);

        // user input event checking logic
        SDL_Event event;
        runEventCheck(&event, &wp, &gb, &sc, &buttons, &stats_window);

        // render the bodies
        body_renderOrbitBodies(renderer, &gb, wp);

        // render the spacecraft
        craft_renderCrafts(renderer, &sc);

        // render stats in main window if enabled
        if (stats_window.is_shown) {
            renderStatsBox(renderer, &gb, &sc, wp, &stats_window);
        }

        // unlock sim vars mutex when done
        pthread_mutex_unlock(&sim_vars_mutex);

        // shows and limits FPS
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

    // free all bodies
    for (int i = 0; i < gb.count; i++) {
        free(gb.names[i]);
    }
    free(gb.names);
    free(gb.mass);
    free(gb.radius);
    free(gb.pixel_radius);
    free(gb.pos_x);
    free(gb.pos_y);
    free(gb.pixel_coordinates_x);
    free(gb.pixel_coordinates_y);
    free(gb.vel_x);
    free(gb.vel_y);
    free(gb.vel);
    free(gb.acc_x);
    free(gb.acc_y);
    free(gb.acc_x_prev);
    free(gb.acc_y_prev);
    free(gb.force_x);
    free(gb.force_y);
    free(gb.kinetic_energy);

    // free all spacecraft
    for (int i = 0; i < sc.count; i++) {
        free(sc.names[i]);
        free(sc.burn_properties[i]);
    }
    free(sc.names);
    free(sc.current_total_mass);
    free(sc.dry_mass);
    free(sc.fuel_mass);
    free(sc.pos_x);
    free(sc.pos_y);
    free(sc.pixel_coordinates_x);
    free(sc.pixel_coordinates_y);
    free(sc.attitude);
    free(sc.vel_x);
    free(sc.vel_y);
    free(sc.vel);
    free(sc.rotational_v);
    free(sc.momentum);
    free(sc.acc_x);
    free(sc.acc_y);
    free(sc.acc_x_prev);
    free(sc.acc_y_prev);
    free(sc.rotational_a);
    free(sc.moment_of_inertia);
    free(sc.grav_force_x);
    free(sc.grav_force_y);
    free(sc.torque);
    free(sc.thrust);
    free(sc.mass_flow_rate);
    free(sc.specific_impulse);
    free(sc.throttle);
    free(sc.nozzle_gimbal_range);
    free(sc.nozzle_velocity);
    free(sc.engine_on);
    free(sc.num_burns);
    free(sc.burn_properties);

    if (g_font) TTF_CloseFont(g_font);
    if (g_font_small) TTF_CloseFont(g_font_small);
    TTF_Quit();
    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();
    return 0;
}