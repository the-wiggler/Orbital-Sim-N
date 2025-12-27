#ifndef TYPES_H
#define TYPES_H

#include <stdbool.h>
#include <stdio.h>
#include <SDL3/SDL.h>
#include <GL/glew.h>
#include <GL/gl.h>

#define PATH_CACHE_LENGTH 100

typedef struct {
    float x, y, z;
} coord_t;

typedef struct {
    int screen_width, screen_height;
    double time_step;
    float window_size_x, window_size_y;

    // 3D camera
    coord_t camera_pos;    // camera position in world space (defined by a unit vector, whereas the magnitude is changed by the viewport zoom)
    float zoom;             // zoom level

    float font_size;
    bool window_open;
    bool data_logging_enabled;
    bool sim_running;
    double sim_time;
    SDL_WindowID main_window_ID;

    double meters_per_pixel;

    bool is_dragging;
    float drag_last_x, drag_last_y;

    bool reset_sim;
    bool is_zooming;
    bool is_zooming_out;
    bool is_zooming_in;

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
    double* pos_z;
    double* vel_x;
    double* vel_y;
    double* vel_z;
    double* vel;
    double* acc_x;
    double* acc_y;
    double* acc_z;
    double* acc_x_prev;
    double* acc_y_prev;
    double* acc_z_prev;
    double* force_x;
    double* force_y;
    double* force_z;
    double* kinetic_energy;
} body_properties_t;

typedef struct {
    bool tangent;   // the burn axis heading beings tangent to the orbit
    bool normal;    // the burn axis heading begins normal to the orbit
    bool absolute;  // the burn axis heading is relative to the vertical axis of space
} relative_burn_target_t;

typedef struct {
    double burn_start_time;
    double burn_end_time;
    double throttle;
    double burn_heading;
    int burn_target_id;
    relative_burn_target_t relative_burn_target; // the axis of rotation the burn heading will be measured from
} burn_properties_t;

typedef struct {
    int count;

    char** names;

    double* current_total_mass;
    double* dry_mass;
    double* fuel_mass;

    double* pos_x;
    double* pos_y;
    double* pos_z;
    double* attitude;

    double* vel_x;
    double* vel_y;
    double* vel_z;
    double* vel;
    double* rotational_v;

    double* momentum;

    double* acc_x;
    double* acc_y;
    double* acc_z;
    double* acc_x_prev;
    double* acc_y_prev;
    double* acc_z_prev;
    double* rotational_a;

    double* moment_of_inertia;

    double* grav_force_x;
    double* grav_force_y;
    double* grav_force_z;
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

// container for all the sim elements
typedef struct {
    body_properties_t gb; // global bodies
    spacecraft_properties_t gs; // global spacecraft
    window_params_t wp; // window properties
    double system_kinetic_energy, system_potential_energy; // total energies of the whole system (reset each iteration)
} sim_properties_t;

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
} button_t;

typedef struct {
    button_t sc_button;
    button_t csv_load_button;
    button_t craft_view_button;
    button_t show_stats_button;
} button_storage_t;


typedef struct {
    FILE* body_pos_FILE;
    FILE* global_data_FILE;
} binary_filenames_t;

typedef struct {
    double timestamp;
    int body_index;
    double pos_data_x, pos_data_y;
    double vel_data_x, vel_data_y;
    double acc_data_x, acc_data_y;
    double force_data_x, force_data_y;
} global_data_t;

typedef struct {
    SDL_Window* window;
    SDL_GLContext glContext;
} SDL_GL_init_t;

typedef struct {
    GLuint VAO;
    GLuint VBO;
} VBO_t;

typedef struct {
    float m[16]; // 4x4 matrix
} mat4;

#endif
