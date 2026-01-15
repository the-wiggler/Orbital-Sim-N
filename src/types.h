#ifndef TYPES_H
#define TYPES_H

#include <stdbool.h>
#include <stdio.h>
#include <SDL3/SDL.h>
#include <GL/glew.h>

#ifdef __APPLE__
#include <OpenGL/gl.h>
#else
#include <GL/gl.h>
#endif

typedef struct {
    float x, y, z;
} vec3_f;

typedef struct {
    double x, y, z;
} vec3;

typedef struct {
    float m[16]; // 4x4 matrix
} mat4;

typedef struct {
    double w, x, y, z;
} quaternion_t;

typedef struct {
    int screen_width, screen_height;
    double time_step;
    float window_size_x, window_size_y;

    // 3D camera
    vec3_f camera_pos;    // camera position in world space (defined by a unit vector, whereas the magnitude is changed by the viewport zoom)
    float zoom;             // zoom level

    volatile bool window_open;
    bool data_logging_enabled;
    volatile bool sim_running;
    double sim_time;
    SDL_WindowID main_window_ID;

    double meters_per_pixel;

    int planet_model_vertex_count;
    int frame_counter;

    bool is_dragging;
    float drag_last_x, drag_last_y;

    bool reset_sim;
    bool is_zooming;
    bool is_zooming_out;
    bool is_zooming_in;

    // visual stuff
    bool draw_lines_between_bodies;
    bool draw_inclination_height;
    bool draw_planet_path;
    bool draw_craft_path;

} window_params_t;

typedef struct {
    char cmd_text_box[256];
    char log[256];
    int cmd_text_box_length;
    float cmd_pos_x, cmd_pos_y;
    float log_pos_x, log_pos_y;
} console_t;

// celestial body
typedef struct {
    char* name;

    double mass;
    double radius;
    double SOI_radius;
    float pixel_radius;

    vec3 pos;
    vec3 vel;
    double vel_mag;
    vec3 acc;
    vec3 acc_prev;
    vec3 force;

    double kinetic_energy;

    double rotational_v;     // angular velocity (rad/s)
    quaternion_t attitude;   // orientation quaternion
} body_t;

// container for all bodies
typedef struct {
    int count;
    int capacity;
    body_t* bodies;
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

// craft
typedef struct {
    char* name;

    double current_total_mass;
    double dry_mass;
    double fuel_mass;

    vec3 pos;
    vec3 vel;
    double vel_mag;
    vec3 acc;
    vec3 acc_prev;
    vec3 grav_force;

    quaternion_t attitude;
    double rotational_v;
    double rotational_a;
    double momentum;
    double moment_of_inertia;
    double torque;

    double thrust;
    double mass_flow_rate;
    double specific_impulse;
    double throttle;
    double nozzle_gimbal_range;
    double nozzle_velocity;
    bool engine_on;

    int SOI_planet_id;
    int SOI_planet_dist;

    int num_burns;
    burn_properties_t* burn_properties;
} spacecraft_t;

// container for all spacecraft
typedef struct {
    int count;
    int capacity;
    spacecraft_t* spacecraft;
} spacecraft_properties_t;

// container for all the sim elements
typedef struct {
    body_properties_t gb; // global bodies
    spacecraft_properties_t gs; // global spacecraft
    window_params_t wp; // window properties
    console_t console; // in-window console
    double system_kinetic_energy, system_potential_energy; // total energies of the whole system (reset each iteration)
} sim_properties_t;

typedef struct {
    int frame_counter;
    bool is_shown;
    double initial_total_energy;
    bool measured_initial_energy;
    double previous_total_energy;

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
    VBO_t vbo;
    float* vertices;
    size_t capacity;  // max number of lines
    size_t count;     // current number of lines
} line_batch_t;

typedef struct {
    float* vertices;
    size_t vertex_count;
    size_t data_size;
} sphere_mesh_t;

// text rendering
typedef struct {
    GLuint tex, shader, vao, vbo;
    float* verts;
    int count;
} font_t;

// planet path tracking
typedef struct {
    vec3* positions;
    int* counts;
    int capacity;
    int num_objects;
} object_path_storage_t;

#endif
