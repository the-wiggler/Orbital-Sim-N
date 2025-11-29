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

// orbital bodies
typedef struct {
    int count;
    int capacity;

    char** names;
    double* mass;
    double* radius;
    float* pixel_radius;
    double* pos_x;
    double* pos_y;
    float* pixel_coordinates_x;
    float* pixel_coordinates_y;
    double* vel_x;
    double* vel_y;
    double* vel;
    double* acc_x;
    double* acc_y;
    double* acc_x_prev; // previous acceleration for verlet integration
    double* acc_y_prev;
    double* force_x; // the forces acting on such body
    double* force_y;
    double* kinetic_energy;
} body_properties_t;

typedef struct {
    double burn_start_time; // simulation time to start burn (seconds)
    double burn_end_time;
    double throttle; // how much throttle you're burning with
    double burn_heading; // (rad)
} burn_properties_t;

// spacecraft
typedef struct {
    int count;
    int capacity;

    char** names;

    double* current_total_mass; // tracks the amount of mass in the ship at an instant
    double* dry_mass; // mass of the empty ship
    double* fuel_mass; // mass of the fuel

    double* pos_x; // in-world coordinates relative to the origin (meters)
    double* pos_y;
    float* pixel_coordinates_x; // on screen coordinates relative to the origin (pixels, center of the screen)
    float* pixel_coordinates_y;
    double* attitude; // 2d heading (radians)

    double* vel_x; // component of velocity relative to space
    double* vel_y;
    double* vel; // total velocity relative to space
    double* rotational_v; // rotational velocity (CCW is positive)

    double* momentum;

    double* acc_x; // component of acceleration relative to space
    double* acc_y;
    double* acc_x_prev; // previous acceleration for verlet integration
    double* acc_y_prev;
    double* rotational_a; // rotational acceleration

    double* moment_of_inertia;

    double* grav_force_x; // sum of the force acting on the body due to the planets' gravitational influence
    double* grav_force_y;
    double* torque; // applied torque (CCW is positive)

    double* thrust; // thrust force of the craft (N)
    double* mass_flow_rate;
    double* specific_impulse;
    double* throttle;
    double* nozzle_gimbal_range; // radians
    double* nozzle_velocity;
    bool* engine_on;

    int* num_burns;
    burn_properties_t** burn_properties; // array of pointers to burn properties for each spacecraft
} spacecraft_properties_t;

typedef struct {
    bool is_shown;
    double initial_total_energy;
    bool measured_initial_energy;
} stats_window_t;

typedef struct {
    body_properties_t* gb;
    spacecraft_properties_t* sc;
    window_params_t* wp;
} physics_sim_args;



void readBodyJSON(const char* FILENAME, body_properties_t* gb);
void readSpacecraftJSON(const char* FILENAME, spacecraft_properties_t* sc);

#endif