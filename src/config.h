#ifndef CONFIG_H
#define CONFIG_H

#include <stdbool.h>
#include <SDL3/SDL.h>
#include <SDL3_ttf/SDL_ttf.h>

typedef struct {
    double time_step;
    int window_size_x, window_size_y;
    int screen_origin_x, screen_origin_y;
    double meters_per_pixel;
    int font_size;
    bool window_open;
    bool sim_running;
    double sim_time;
    SDL_WindowID main_window_ID;
} window_params_t;

typedef struct {
    char* name;
    double mass;
    int radius; // this is an approximate value just for visual purposes
    int pixel_radius;
    double r_from_body;
    double pos_x;
    double pos_y;
    int pixel_coordinates_x;
    int pixel_coordinates_y;
    double vel_x;
    double vel_y;
    double vel;
    double acc_x;
    double acc_y;
    double acc_x_prev; // previous acceleration for Verlet integration
    double acc_y_prev;
    double force_x; // the forces acting on such body
    double force_y;
    double kinetic_energy;
} body_properties_t;

typedef struct {
    char* name;
    double pos_x;
    double pos_y;
    int pixel_coordinates_x;
    int pixel_coordinates_y;
    double vel_x;
    double vel_y;
    double vel;
    double acc_x;
    double acc_y;
    double acc_x_prev; // previous acceleration for verlet integration
    double acc_y_prev;
    double force_x; // sum of the force acting on the body due to the planets' gravitational influence
    double force_y;

    double heading;
    double dry_mass;
    double fuel_mass;
    double current_total_mass; // tracks the amount of mass in the ship at an instant
    double mass_flow_rate;
    double thrust;
    double specific_impulse;
    float throttle;
    bool engine_on;

    double burn_start_time;   // simulation time to start burn (seconds)
    double burn_duration;     // how long to burn
    double burn_heading;      // direction to burn - rad
    double burn_throttle;     // throttle setting for the burn - 0-1
} spacecraft_properties_t;

typedef struct {
    SDL_Window* window;
    SDL_Renderer* renderer;
    SDL_WindowID window_ID;
    bool is_shown;
} stats_window_t;

typedef struct {
    body_properties_t* gb;
} body_pos_storage;

#endif