#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <SDL3/SDL.h>
#include <SDL3/SDL_main.h>
#include <SDL3_ttf/SDL_ttf.h>
#include <string.h>
#include <stdbool.h>
#include <time.h>
#include <unistd.h>
#include "calculation_functions.h"
#include "sdl_elements.h"

// NOTE: ALL CALCULATIONS SHOULD BE DONE IN BASE SI UNITS

// universal gravitation constant
const double G = 6.67430E-11;
TTF_Font* g_font = NULL;

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
    wp.speed_multiplier = 1.0; // multiplier for TIME_STEP
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

    speed_control_t speed_control = {
        wp.window_size_x * 0.01,  // x position
        wp.window_size_y * 0.01,  // y position
        wp.window_size_x * 0.15,  // width
        wp.window_size_y * 0.04,  // height
        false                                // is_hovered
    };

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
    // SDL ttf font stuff
    TTF_Init();
    g_font = TTF_OpenFont("CascadiaCode.ttf", wp.font_size);

    addOrbitalBody(&global_bodies, &num_bodies, 1e24, 3e8, 0, 0, -1000);
    addOrbitalBody(&global_bodies, &num_bodies, 5e24, 0, 0, 0, 0);

    ////////////////////////////////////////////////////////
    // simulation loop                                    //
    ////////////////////////////////////////////////////////
    while (wp.window_open) {
        // checks inputs into the window
        SDL_Event event;
        runEventCheck(&event, &speed_control, &wp, &global_bodies, &num_bodies);

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
        for (int i = 0; i < num_bodies; i++) {
            // draw bodies
            SDL_RenderFillCircle(renderer, global_bodies[i].pixel_coordinates_x,
                            global_bodies[i].pixel_coordinates_y, 
                            calculateVisualRadius(global_bodies[i], wp));
        }

        ////////////////////////////////////////////////////
        // UI ELEMENTS                                    //
        ////////////////////////////////////////////////////
        // draw scale reference bar
        drawScaleBar(renderer, wp);

        // draw speed control box
        drawSpeedControl(renderer, &speed_control, wp);

        // draw stats box
        if (global_bodies != NULL) drawStatsBox(renderer, global_bodies, num_bodies, wp.sim_time, wp);

        // help text at the bottom
        SDL_WriteText(renderer, g_font, "Space: pause/resume | Left Click: add body | R: Reset", wp.window_size_x * 0.4, wp.window_size_y - wp.window_size_x * 0.02 - wp.font_size, white_text);
        // pause indicator
        if (wp.sim_running) SDL_WriteText(renderer, g_font, "Sim Running...", wp.window_size_x * 0.75, wp.window_size_y * 0.015, white_text);
        else SDL_WriteText(renderer, g_font, "Sim Paused", wp.window_size_x * 0.75, wp.window_size_y * 0.015, white_text);

        // present the renderer to the screen
        SDL_RenderPresent(renderer);
    }

    // clean up
    free(global_bodies);
    global_bodies = NULL;
    num_bodies = 0;
    if (g_font) TTF_CloseFont(g_font);
    TTF_Quit();
    SDL_DestroyWindow(window);
    SDL_Quit();
    return 0;
}