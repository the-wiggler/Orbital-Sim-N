#include "config.h"
#include "stats_window.h"
#ifdef _WIN32
    #include <windows.h>
#else
    #include <unistd.h>
#endif
#include <threads.h>
#include <stdlib.h>
#include <math.h>
#include <SDL3/SDL.h>
#include <SDL3/SDL_main.h>
#include <SDL3_ttf/SDL_ttf.h>
#include <stdbool.h>
#include "sim_calculations.h"
#include "sdl_elements.h"

// NOTE: ALL CALCULATIONS SHOULD BE DONE IN BASE SI UNITS

// universal gravitation constant
char* FILENAME = "planet_data.csv";
char* SPACECRAFT_FILENAME = "spacecraft_data.csv";
const double G = 6.67430E-11;
TTF_Font* g_font = NULL;
TTF_Font* g_font_small = NULL;

////////////////////////////////////////////////////////////////////////////////////////////////////
// MAIN
////////////////////////////////////////////////////////////////////////////////////////////////////
int main(int argc, char* argv[]) {
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
    SDL_Window* window = SDL_CreateWindow("Orbit Simulation", wp.window_size_x, wp.window_size_y, SDL_WINDOW_RESIZABLE);
    wp.main_window_ID = SDL_GetWindowID(window);
    // create an SDL renderer and clear the window to create a blank canvas
    SDL_Renderer *renderer = SDL_CreateRenderer(window, NULL);
    SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
    SDL_RenderClear(renderer);

    // SDL ttf font stuff
    TTF_Init();
    g_font = TTF_OpenFont("CascadiaCode.ttf", wp.font_size);
    g_font_small = TTF_OpenFont("CascadiaCode.ttf", (float)wp.window_size_x / 90);


    // toggleable stats window
    stats_window_t stats_window = {0};
    SDL_Color white_text = {255, 255, 255, 255};

    ////////////////////////////////////////
    // SIM VARS                           //
    ////////////////////////////////////////

    // initialize simulation objects
    int num_bodies = 0;
    body_properties_t* gb = NULL;
    int num_craft = 0;
    spacecraft_properties_t* sc = NULL;

    // initialize planet position storage variable
    // this variable should be used as an array of gb arrays (im sorry)
    // ***CAREFUL!!! this variable is shared between the window and sim threads
    body_pos_storage* b_pos_data = NULL;
    mtx_t b_pod_data_mutex;

    ////////////////////////////////////////////////////////
    // simulation loop                                    //
    ////////////////////////////////////////////////////////
    while (wp.window_open) {
        // user input event checking logic
        SDL_Event event;
        runEventCheck(&event, &wp, &gb, &num_bodies, &sc, &num_craft, &buttons, &dialog, &stats_window);

        // clears previous frame from the screen
        SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
        SDL_RenderClear(renderer);

        // draw scale reference bar
        drawScaleBar(renderer, wp);
        
        ////////////////////////////////////////////////////
        // START OF SIMULATION LOGIC                      //
        ////////////////////////////////////////////////////
        // IMPORTANT -- DOES ALL OF THE BODY CALCULATIONS:
        runCalculations(&gb, &sc, &wp, num_bodies, num_craft);

        // render the bodies
        body_renderOrbitBodies(renderer, gb, num_bodies, wp);

        // render the spacecraft
        craft_renderCrafts(renderer, sc, num_craft, wp);

        ////////////////////////////////////////////////////
        // UI ELEMENTS                                    //
        ////////////////////////////////////////////////////

        // draw speed control button
        renderUIButtons(renderer, &buttons, &wp);

        // help text at the bottom
        SDL_WriteText(renderer, g_font, "Space: pause/resume | R: Reset", wp.window_size_x * 0.4, wp.window_size_y - wp.window_size_x * 0.02 - wp.font_size, white_text);

        // render text input dialog if active
        renderBodyTextInputDialog(renderer, &dialog, wp);

        // draw time indicator text
        renderTimeIndicators(renderer, wp);

        // render the stats window if active
        if (stats_window.is_shown) StatsWindow_render(&stats_window, 60, 0, 0, gb, num_bodies, wp);

        // present the renderer to the screen
        SDL_RenderPresent(renderer);
    }
    ////////////////////////////////////////////////////////
    // end of simulation loop                             //
    ////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////
    // CLEAN UP                                       //
    ////////////////////////////////////////////////////
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

    StatsWindow_destroy(&stats_window);
    if (g_font) TTF_CloseFont(g_font);
    if (g_font_small) TTF_CloseFont(g_font_small);
    TTF_Quit();
    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();
    return 0;
}