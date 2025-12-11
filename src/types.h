#ifndef TYPES_H
#define TYPES_H

#include <stdbool.h>
#include <SDL3/SDL.h>

#define PATH_CACHE_LENGTH 100

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

    bool reset_sim;
    bool is_zooming;

    bool craft_view_shown;
    bool main_view_shown;
} window_params_t;

typedef struct {
    int count;

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
    double* acc_x_prev;
    double* acc_y_prev;
    double* force_x;
    double* force_y;
    double* kinetic_energy;

    int cache_counter;
    int cache_valid_count;
    float** cached_body_coords_x;
    float** cached_body_coords_y;
} body_properties_t;

typedef struct {
    double burn_start_time;
    double burn_end_time;
    double throttle;
    double burn_heading;
} burn_properties_t;

typedef struct {
    int count;

    char** names;

    double* current_total_mass;
    double* dry_mass;
    double* fuel_mass;

    double* pos_x;
    double* pos_y;
    float* pixel_coordinates_x;
    float* pixel_coordinates_y;
    double* attitude;

    double* vel_x;
    double* vel_y;
    double* vel;
    double* rotational_v;

    double* momentum;

    double* acc_x;
    double* acc_y;
    double* acc_x_prev;
    double* acc_y_prev;
    double* rotational_a;

    double* moment_of_inertia;

    double* grav_force_x;
    double* grav_force_y;
    double* torque;

    double* thrust;
    double* mass_flow_rate;
    double* specific_impulse;
    double* throttle;
    double* nozzle_gimbal_range;
    double* nozzle_velocity;
    bool* engine_on;

    int* num_burns;
    burn_properties_t** burn_properties;
} spacecraft_properties_t;

typedef struct {
    int frame_counter;
    bool is_shown;
    double initial_total_energy;
    bool measured_initial_energy;
    double previous_total_energy;

    char vel_text[10][64];
    char total_energy_text[64];
    char delta_energy_text[64];
    char error_text[64];
    SDL_Color error_color;
    int cached_body_count;
} stats_window_t;

typedef struct {
    float x, y, width, height;
    bool is_hovered;
    SDL_Color normal_color;
    SDL_Color hover_color;
    SDL_Texture* normal_texture;
    SDL_Texture* hover_texture;
    bool textures_valid;
} button_t;

typedef struct {
    button_t sc_button;
    button_t csv_load_button;
    button_t craft_view_button;
    button_t show_stats_button;
} button_storage_t;

typedef struct {
    body_properties_t* gb;
    spacecraft_properties_t* sc;
    window_params_t* wp;
} physics_sim_args;

typedef struct {
    body_properties_t* gb;
    spacecraft_properties_t* sc;
    window_params_t* wp;
} cleanup_args;

#endif
