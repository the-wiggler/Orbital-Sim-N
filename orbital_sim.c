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

// NOTE: ALL CALCULATIONS SHOULD BE DONE IN BASE SI UNITS

////////////////////////////////////////////////////////////////////////////////////////////////////
// CONSTANTS AND GLOBAL VARIABLES
////////////////////////////////////////////////////////////////////////////////////////////////////
///////////// simulation constants:
// universal gravitation constant
const double G = 6.67430E-11;

double TIME_STEP                = 1; // amount of time in seconds each step in the simulation should be
                                        // i.e. each new updated position shown is after x seconds
                                        // this value can be changed to adjust the speed/accuracy of the simulation
speed_control_t speed_control = {10, 10, 150, 40, false};
double speed_multiplier = 1.0; // multiplier for TIME_STEP

// SDL window sizing numbers
const int WINDOW_SIZE_X         = 1000;
const int WINDOW_SIZE_Y         = 1000;
const int ORIGIN_X              = WINDOW_SIZE_X / 2;
const int ORIGIN_Y              = WINDOW_SIZE_Y / 2;
double meters_per_pixel   = 100000; // 1m in space will equal x number of pixels on screen
const int FONT_SIZE             = WINDOW_SIZE_Y / (WINDOW_SIZE_X * 0.05);

TTF_Font* g_font = NULL;

// holds the information on the bodies
int num_bodies = 0;
body_properties_t* global_bodies = NULL;

////////////////////////////////////////////////////////////////////////////////////////////////////
// MAIN
////////////////////////////////////////////////////////////////////////////////////////////////////
int main(int argc, char* argv[]) {
    // initialize SDL3
    SDL_Init(SDL_INIT_VIDEO);
    // create an SDL window
    SDL_Window* window = SDL_CreateWindow("Orbit Simulation", WINDOW_SIZE_X, WINDOW_SIZE_Y, SDL_WINDOW_RESIZABLE);
    // create an SDL renderer and clear the window to create a blank canvas
    SDL_Renderer *renderer = SDL_CreateRenderer(window, NULL);
    SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
    SDL_RenderClear(renderer);
    // SDL ttf font stuff
    TTF_Init();
    g_font = TTF_OpenFont("CascadiaCode.ttf", FONT_SIZE);

    // set to false to stop the sim program
    bool sim_running = true;

    ////////////////////////////////////////////////////////
    // simulation loop                                    //
    ////////////////////////////////////////////////////////

    // Earth at origin
    addOrbitalBody(5.972e24, 0, 0, 0, 0);

    // Moon orbiting Earth
    addOrbitalBody(7.342e22, 3.844e8, 0, 0, 1022);

    addOrbitalBody(7.342e22, -3.844e8, 0, 0, -1022);

    while (sim_running) {
        // checks inputs into the window
        SDL_Event event;
        runEventCheck(&event, &sim_running, &speed_control, &TIME_STEP, &meters_per_pixel);

        // runs a check to see what the zoom value should be to automatically asjust

        // clears previous frame from the screen
        SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
        SDL_RenderClear(renderer);

        // color for drawing bodies
        SDL_SetRenderDrawColor(renderer, 255, 255, 255, 255);
        
        ////////////////////////////////////////////////////
        // START OF SIMULATION LOGIC                      //
        ////////////////////////////////////////////////////

        // calculate forces between all body pairs
        if (global_bodies != NULL) {
            for (int i = 0; i < num_bodies; i++) {
                global_bodies[i].force_x = 0;
                global_bodies[i].force_y = 0;
                for (int j = 0; j < num_bodies; j++) {
                    if (i != j) {
                        calculateForce(&global_bodies[i], global_bodies[j]);
                    }
                }
            }
            
            // update the motion for each body and draw
            for (int i = 0; i < num_bodies; i++) {
                // updates the kinematic properties of each body (velocity, accelertion, position, etc)
                updateMotion(&global_bodies[i], TIME_STEP);
                // transform real-space coordinate to pixel coordinates on screen (scaling)
                transformCoordinates(&global_bodies[i]);
                // draw bodies
                SDL_RenderFillCircle(renderer, global_bodies[i].pixel_coordinates_x,
                                global_bodies[i].pixel_coordinates_y, 
                                calculateVisualRadius(global_bodies[i]));
            }
        }


        // draw scale reference bar
        drawScaleBar(renderer, meters_per_pixel, WINDOW_SIZE_X, WINDOW_SIZE_Y);

        // draw speed control box
        drawSpeedControl(renderer, &speed_control, TIME_STEP);

        // draw stats box
        if (global_bodies != NULL) drawStatsBox(renderer, global_bodies, num_bodies);

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