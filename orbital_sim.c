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
speed_control_t speed_control   = {10, 10, 150, 40, false};
double speed_multiplier = 1.0; // multiplier for TIME_STEP

// SDL window sizing numbers
const int WINDOW_SIZE_X         = 1500;
const int WINDOW_SIZE_Y         = 1500;
const int ORIGIN_X              = WINDOW_SIZE_X / 2;
const int ORIGIN_Y              = WINDOW_SIZE_Y / 2;
double meters_per_pixel         = 100000; // 1m in space will equal x number of pixels on screen
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
    bool window_open = true;
    bool sim_running = true;
    double sim_time  = 0;

    ////////////////////////////////////////////////////////
    // simulation loop                                    //
    ////////////////////////////////////////////////////////

    double jupiterMass = 1.898e27;   // Jupiter
    double ioMass = 8.93e22;         // Io (most volcanically active body)
    double ioRadius = 4.217e8;       // Io's orbital radius
    double ioVel = 17334.0;  // m/s
    double europaMass = 4.8e22;      // Europa
    double europaRadius = 6.709e8;   // Europa's orbital radius
    double europaVel = 13740.0;  // m/s
    
    addOrbitalBody(jupiterMass, 0.0, 0.0, 0.0, 0.0);
    addOrbitalBody(ioMass, ioRadius, 0.0, 0.0, ioVel);
    addOrbitalBody(europaMass, 0.0, europaRadius, -europaVel, 0.0);
    addOrbitalBody(1e25, 1e8, 1e8, 27334, 0);

    while (window_open) {
        // checks inputs into the window
        SDL_Event event;
        runEventCheck(&event, &window_open, &speed_control, &TIME_STEP, &meters_per_pixel);

        // clears previous frame from the screen
        SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
        SDL_RenderClear(renderer);

        // color for drawing bodies
        SDL_SetRenderDrawColor(renderer, 255, 255, 255, 255);
        
        ////////////////////////////////////////////////////
        // START OF SIMULATION LOGIC                      //
        ////////////////////////////////////////////////////
        if (sim_running) {
            if (sim_time > 10000) sim_running = false;
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
            sim_time += TIME_STEP;
        }

        // draw scale reference bar
        drawScaleBar(renderer, meters_per_pixel, WINDOW_SIZE_X, WINDOW_SIZE_Y);

        // draw speed control box
        drawSpeedControl(renderer, &speed_control, TIME_STEP);

        // draw stats box
        if (global_bodies != NULL) drawStatsBox(renderer, global_bodies, num_bodies, sim_time);

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