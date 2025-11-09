#include <stdio.h>
#include <math.h>
#include <SDL3/SDL.h>
#include <string.h>
#include <stdbool.h>
#include <time.h>
#include <unistd.h>

// NOTE: ALL CALCULATIONS SHOULD BE DONE IN BASE SI UNITS

////////////////////////////////////////////////////////////////////////////////////////////////////
// CONSTANTS
////////////////////////////////////////////////////////////////////////////////////////////////////
///////////// simulation constants:
// universal gravitation constant
const double G = 6.67430E-11;
// fixed time step
const double TIME_STEP = 1000; // amount of time in seconds each step in the simulation should be
                                // i.e. each new updated position shown is after x seconds

// SDL window sizing numbers
const int WINDOW_SIZE_X = 1000;
const int WINDOW_SIZE_Y = 1000;
const int ORIGIN_X = WINDOW_SIZE_X / 2;
const int ORIGIN_Y = WINDOW_SIZE_Y / 2;
const double METERS_PER_PIXEL = 1000; // 1m in space will equal x number of pixels on screen

////////////////////////////////////////////////////////////////////////////////////////////////////
// STRUCTS
////////////////////////////////////////////////////////////////////////////////////////////////////
typedef struct {
    double mass;
    double radius;
    double pos_x;
    double pos_y;
    double vel_x;
    double vel_y;
    double acc_x;
    double acc_y;
    double force_x;
    double force_y;
    int pixel_coordinates_x;
    int pixel_coordinates_y;
} body_properties_t;

////////////////////////////////////////////////////////////////////////////////////////////////////
// UNIVERSAL VARIABLES
////////////////////////////////////////////////////////////////////////////////////////////////////
body_properties_t body1 = {0};
body_properties_t body2 = {0};

////////////////////////////////////////////////////////////////////////////////////////////////////
// SIMULATION CALCULATION FUNCTIONS
////////////////////////////////////////////////////////////////////////////////////////////////////
// at the end of each sim loop, this function should be run to calculate the changes in
// the force values based on other parameters. for example, using F to find a based on m.
void calculateForces(body_properties_t *b, body_properties_t b2) {
    // calculate the distance between the two bodies
    double delta_pos_x = b2.pos_x - b->pos_x;
    double delta_pos_y = b2.pos_y - b->pos_y;
    double r = sqrt(delta_pos_x * delta_pos_x + delta_pos_y * delta_pos_y);

    // calculate the force that b2 applies on b due to gravitation (F = (GMm) / r)
    double total_force = (G * b->mass * b2.mass) / (r * r);
    b->force_x = total_force * (delta_pos_x / r);
    b->force_y = total_force * (delta_pos_y / r);

    // calculate the acceleration from the force on the object
    b->acc_x = b->force_x / b->mass;
    b->acc_y = b->force_y / b->mass;
}

// this calculates the changes of velocity and position based on the force values from before
void updateMotion(body_properties_t *b, double dt) {
    // update the velocity
    b->vel_x = b->vel_x + b->acc_x * dt;
    b->vel_y = b->vel_y + b->acc_y * dt;

    // update the position using the new velocity
    b->pos_x = b->pos_x + b->vel_x * dt;
    b->pos_y = b->pos_y + b->vel_y * dt;
}

// transforms spacial coordinates (for example, in meters) to pixel coordinates
void transformCoordinates(body_properties_t *b) {
    b->pixel_coordinates_x = ORIGIN_X + (int)(b->pos_x / METERS_PER_PIXEL);
    b->pixel_coordinates_y = ORIGIN_Y + (int)(b->pos_y / METERS_PER_PIXEL);
}

// draw a circle in SDL
void SDL_RenderFillCircle(SDL_Renderer* renderer, int centerX, int centerY, int radius) {
    int x = radius;
    int y = 0;
    int radiusError = 1 - x;

    while (x >= y) {
        SDL_RenderLine(renderer, centerX - x, centerY + y, centerX + x, centerY + y);
        SDL_RenderLine(renderer, centerX - x, centerY - y, centerX + x, centerY - y);
        SDL_RenderLine(renderer, centerX - y, centerY + x, centerX + y, centerY + x);
        SDL_RenderLine(renderer, centerX - y, centerY - x, centerX + y, centerY - x);

        y++;
        if (radiusError < 0) {
            radiusError += 2 * y + 1;
        } else {
            x--;
            radiusError += 2 * (y - x + 1);
        }
    }
}

// draw a distance scale bar on the sreen
void drawScaleBar(SDL_Renderer* renderer, double meters_per_pixel, int window_width, int window_height) {
    const int BAR_HEIGHT = 3;
    const int MARGIN = 20;
    const double REFERENCE_DISTANCE_KM = 10.0;
    const double REFERENCE_DISTANCE_M  = REFERENCE_DISTANCE_KM * 1000;

    int bar_width_pixels = (int)(REFERENCE_DISTANCE_M / meters_per_pixel);

    // position of the scale bar
    int bar_x = MARGIN;
    int bar_y = window_height - MARGIN - BAR_HEIGHT;

    // draw the bar
    SDL_SetRenderDrawColor(renderer, 255, 255, 255, 255);
    SDL_FRect bar_rect = {bar_x, bar_y, bar_width_pixels, BAR_HEIGHT};
    SDL_RenderFillRect(renderer, &bar_rect);
    SDL_FRect left_cap = {bar_x, bar_y - 3, 2, BAR_HEIGHT + 6};
    SDL_FRect right_cap = {bar_x + bar_width_pixels - 2, bar_y - 3, 2, BAR_HEIGHT + 6};
    SDL_RenderFillRect(renderer, &left_cap);
    SDL_RenderFillRect(renderer, &right_cap);
}

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

    // set to false to stop the sim program
    bool sim_running = true;

    // initial properties of the bodies
    body1.mass = 500000000000.0f;
    body2.mass = 500000000000.0f;
    body1.pos_y = 25000.0f;
    body2.pos_y = -25000.0f;
    body1.vel_x = 0.01f;
    body2.vel_x = -0.01f;

    ////////////////////////////////////////////////////////
    // simulation loop                                    //
    ////////////////////////////////////////////////////////
    while (sim_running) {
        // checks if I pressed the button to close the program
        SDL_Event event;
        while (SDL_PollEvent(&event)) {
            if (event.type == SDL_EVENT_QUIT) {
                sim_running = false;
            }
        }

        // clears previous frame from the screen
        SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
        SDL_RenderClear(renderer);

        // calculate the forces that each body applies on the other
        calculateForces(&body1, body2); // finds the force body 2 applies on body 1
        calculateForces(&body2, body1); // finds the force body 1 applies on body 2

        // update the velocity and position of both bodies based on the calculated force
        updateMotion(&body1, TIME_STEP);
        updateMotion(&body2, TIME_STEP);

        // map position data to pixel data on the screen
        transformCoordinates(&body1);
        transformCoordinates(&body2);

        // set the renderer color and draw the bodies as circles to the renderer
        SDL_SetRenderDrawColor(renderer, 91, 161, 230, 255);
        SDL_RenderFillCircle(renderer, body1.pixel_coordinates_x, body1.pixel_coordinates_y, 4);
        SDL_SetRenderDrawColor(renderer, 230, 166, 91, 255);
        SDL_RenderFillCircle(renderer, body2.pixel_coordinates_x, body2.pixel_coordinates_y, 4);

        // draw scale reference bar
        drawScaleBar(renderer, METERS_PER_PIXEL, WINDOW_SIZE_X, WINDOW_SIZE_Y);

        // present the renderer to the screen
        SDL_RenderPresent(renderer);
    }

    // clean up
    SDL_DestroyWindow(window);
    SDL_Quit();
    return 0;
}