#ifndef CONFIG_H
#define CONFIG_H

#include <stdbool.h>
#include <SDL3/SDL.h>

typedef struct {
    int screen_width, screen_height;
    double time_step;
    float window_size_x, window_size_y;
    float screen_origin_x, screen_origin_y;
    double meters_per_pixel;
    float font_size;
    bool window_open;
    bool sim_running;
    double sim_time;
    SDL_WindowID main_window_ID;

    bool is_dragging;
    float drag_start_x;
    float drag_start_y;
    float drag_origin_x;
    float drag_origin_y;
} window_params_t;

typedef struct {
    char* name;
    double mass;
    float radius; // this is an approximate value just for visual purposes
    float pixel_radius;
    double pos_x;
    double pos_y;
    float pixel_coordinates_x;
    float pixel_coordinates_y;
    double vel_x;
    double vel_y;
    double vel;
    double acc_x;
    double acc_y;
    double acc_x_prev; // previous acceleration for verlet integration
    double acc_y_prev;
    double force_x; // the forces acting on such body
    double force_y;
    double kinetic_energy;
} body_properties_t;

typedef struct {
    char* name;
    double pos_x;
    double pos_y;
    float pixel_coordinates_x;
    float pixel_coordinates_y;
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
    double throttle;
    bool engine_on;

    double burn_start_time;   // simulation time to start burn (seconds)
    double burn_duration;     // how long to burn
    double burn_heading;      // direction to burn - rad
    double burn_throttle;     // throttle setting for the burn - 0-1
} spacecraft_properties_t;

typedef struct {
    bool is_shown;
    double initial_total_energy;
    bool measured_initial_energy;
} stats_window_t;

typedef struct {
    body_properties_t** gb;
    spacecraft_properties_t** sc;
    window_params_t* wp;
    int* num_bodies;
    int* num_craft;
} physics_sim_args;



void readBodyJSON(const char* FILENAME, body_properties_t** gb, int* num_bodies);
void readSpacecraftJSON(const char* FILENAME, spacecraft_properties_t** sc, int* num_craft);

#endif