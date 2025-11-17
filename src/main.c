#include "config.h"
#include "stats_window.h"
#ifdef _WIN32
    #include <windows.h>
#else
    #include <unistd.h>
#endif
#include <stdlib.h>
#include <math.h>
#include <SDL3/SDL.h>
#include <SDL3/SDL_main.h>
#include <SDL3_ttf/SDL_ttf.h>
#include <stdbool.h>
#include "calculation_functions.h"
#include "sdl_elements.h"

// NOTE: ALL CALCULATIONS SHOULD BE DONE IN BASE SI UNITS

// universal gravitation constant
char* FILENAME = "planet_data.csv";
const double G = 6.67430E-11;
TTF_Font* g_font = NULL;
TTF_Font* g_font_small = NULL;

////////////////////////////////////////////////////////////////////////////////////////////////////
// MAIN
////////////////////////////////////////////////////////////////////////////////////////////////////
int main(int argc, char* argv[]) {
    ////////////////////////////////////////////////////////
    // sim variables                                      //
    ////////////////////////////////////////////////////////
    window_params_t wp = {0};

    wp.time_step = 1; // amount of time in seconds each step in the simulation should be
                                // i.e. each new updated position shown is after x seconds
                                // this value can be changed to adjust the speed/accuracy of the simulation
    // SDL window sizing numbers
    wp.window_size_x = 1000;
    wp.window_size_y = 1000;
    wp.screen_origin_x = wp.window_size_x / 2;
    wp.screen_origin_y = wp.window_size_y / 2;
    wp.meters_per_pixel = 100000;
    wp.font_size = wp.window_size_x / 50;

    // set to false to stop the sim program
    wp.window_open = true;
    wp.sim_running = true; // set to false to pause simulation, set to true to resume
    wp.sim_time  = 0; // tracks the passed time in simulation

    button_storage_t buttons;
    initButtons(&buttons, wp);

    // initialize text input dialog
    text_input_dialog_t dialog = {0};
    dialog.active = false;
    dialog.state = INPUT_NONE;

    SDL_Color white_text = {255, 255, 255, 255};

    // holds the information on the bodies
    int num_bodies;
    body_properties_t* global_bodies = NULL;
    if (global_bodies == NULL) num_bodies = 0;

    // initialize SDL3
    SDL_Init(SDL_INIT_VIDEO);
    // create an SDL window
    SDL_Window* window = SDL_CreateWindow("Orbit Simulation", wp.window_size_x, wp.window_size_y, SDL_WINDOW_RESIZABLE);
    // create an SDL renderer and clear the window to create a blank canvas
    SDL_Renderer *renderer = SDL_CreateRenderer(window, NULL);
    SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
    SDL_RenderClear(renderer);

    // toggleable stats window
    stats_window_t stats_window = {0};

    // SDL ttf font stuff
    TTF_Init();
    g_font = TTF_OpenFont("CascadiaCode.ttf", wp.font_size);
    g_font_small = TTF_OpenFont("CascadiaCode.ttf", (float)wp.window_size_x / 90);

    ////////////////////////////////////////////////////////
    // simulation loop                                    //
    ////////////////////////////////////////////////////////
    while (wp.window_open) {
        // checks inputs into the window
        SDL_Event event;
        runEventCheck(&event, &wp, &global_bodies, &num_bodies, &buttons, &dialog, &stats_window);

        // clears previous frame from the screen
        SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
        SDL_RenderClear(renderer);

        // color for drawing bodies
        SDL_SetRenderDrawColor(renderer, 255, 255, 255, 255);
        
        ////////////////////////////////////////////////////
        // START OF SIMULATION LOGIC                      //
        ////////////////////////////////////////////////////
        // IMPORTANT -- DOES ALL OF THE BODY CALCULATIONS:
        runCalculations(&global_bodies, &wp, num_bodies);

        // render the bodies
        renderOrbitBodies(renderer, global_bodies, num_bodies, wp);

        ////////////////////////////////////////////////////
        // UI ELEMENTS                                    //
        ////////////////////////////////////////////////////
        // draw scale reference bar
        drawScaleBar(renderer, wp);

        // draw speed control button
        renderUIButtons(renderer, &buttons, &wp);

        // help text at the bottom
        SDL_WriteText(renderer, g_font, "Space: pause/resume | R: Reset", wp.window_size_x * 0.4, wp.window_size_y - wp.window_size_x * 0.02 - wp.font_size, white_text);

        // render text input dialog if active
        renderBodyTextInputDialog(renderer, &dialog, wp);

        // render the stats window if active
        if (stats_window.is_shown) {
            StatsWindow_render(&stats_window, 60, 0, 0, global_bodies, num_bodies, wp);
        }

        // present the renderer to the screen
        SDL_RenderPresent(renderer);
    }

    // clean up
    free(global_bodies);
    global_bodies = NULL;
    num_bodies = 0;
    StatsWindow_destroy(&stats_window);
    if (g_font) TTF_CloseFont(g_font);
    if (g_font_small) TTF_CloseFont(g_font_small);
    TTF_Quit();
    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();
    return 0;
}