#ifndef TYPES_H
#define TYPES_H

#include <stdbool.h>
#include <stdio.h>
#ifdef GUI_ENABLED
#include <SDL3/SDL.h>
#ifdef __EMSCRIPTEN__
    #include <GLES3/gl3.h>
#else
    #include <GL/glew.h>
#endif
#ifdef __APPLE__
#include <OpenGL/gl.h>
#else
#include <GL/gl.h>
#endif
#endif
#include "globals.h"

typedef struct {
    float x, y;
} vec2_f;

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
    char cmd_text_box[256];
    char log[256];
    int cmd_text_box_length;
    float cmd_pos_x, cmd_pos_y;
    float log_pos_x, log_pos_y;
} console_t;

typedef struct {
    int SOI_planet_id;
    int closest_planet_id;
    double closest_r_squared;

    double apoapsis, periapsis;

    double semi_major_axis; // m
    double eccentricity;
    double inclination; // rad
    double ascending_node; // rad
    double arg_periapsis; // rad
    double true_anomaly; // rad
    double specific_E; // specific energy
} orbital_elements_t;

// celestial body
typedef struct {
    char name[MAX_NAME_LENGTH];

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

    orbital_elements_t oe;

    double rotational_v;     // angular velocity (rad/s)
    vec3 spin_axis;          // world-space rotation axis (cached)
    quaternion_t delta_rotation_per_step;  // pre-computed rotation quaternion per timestep
    quaternion_t attitude;   // orientation quaternion
} body_t;

// container for all bodies
typedef struct {
    int count;
    double total_system_energy;
    body_t bodies[MAX_PLANETS];
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
    char name[MAX_NAME_LENGTH];

    double current_total_mass; // kg
    double dry_mass; // kg
    double fuel_mass; // kg


    vec3 pos; // meters
    vec3 vel; // m/s
    double vel_mag;
    vec3 acc; // m/s^2
    vec3 acc_prev;
    vec3 grav_force; // N

    quaternion_t attitude; // oh jeez
    double rotational_v; // rad/s
    double rotational_a; // rad/s^2
    double momentum;
    double moment_of_inertia;
    double torque; // N m

    double thrust; // N
    double mass_flow_rate; // kg/s
    double specific_impulse;
    double throttle;
    double nozzle_gimbal_range;
    double nozzle_velocity;
    bool engine_on;

    orbital_elements_t oe;

    int num_burns;
    burn_properties_t burn_properties[MAX_BURNS_PER_SPACECRAFT];
} spacecraft_t;

// container for all spacecraft
typedef struct {
    int count;
    spacecraft_t spacecraft[MAX_SPACECRAFT];
} spacecraft_properties_t;

#ifdef GUI_ENABLED
typedef struct {
    int screen_width, screen_height;
    double time_step;
    float window_size_x, window_size_y;

    // 3D camera
    vec3_f camera_pos;    // camera position in world space (defined by a unit vector, whereas the magnitude is changed by the viewport zoom)
    float zoom;           // zoom level
    float cam_pitch, cam_yaw;
    vec3_f cam_target; // should be the world coordinates

    volatile bool window_open;
    bool data_logging_enabled;
    volatile bool sim_running;
    double sim_time;
    SDL_WindowID main_window_ID;

    double meters_per_pixel;

    int planet_model_vertex_count;
    int frame_counter;

    bool is_dragging_orbit_view;
    vec2_f drag_last;
    bool is_dragging_translation_view; // activating this means you want to edit camera_pos

    bool reset_sim;
    bool is_zooming;
    bool is_zooming_out;
    bool is_zooming_in;

    // visual stuff
    bool draw_lines_between_bodies;
    bool draw_inclination_height;
    bool draw_planet_path;
    bool draw_craft_path;
    bool draw_planet_SOI;

} window_params_t;
#else
typedef struct {
    double time_step;

    volatile bool window_open;
    bool data_logging_enabled;
    volatile bool sim_running;
    double sim_time;

    bool reset_sim;

    double meters_per_pixel;

    int frame_counter;
} window_params_t;
#endif

// container for all the sim elements
typedef struct {
    body_properties_t global_bodies; // global bodies
    spacecraft_properties_t global_spacecraft; // global spacecraft
    window_params_t window_params; // window properties
    console_t console; // in-window console
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

#ifdef GUI_ENABLED
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
    float vertices[MAX_LINE_BATCH * 12];
    size_t capacity;  // max number of lines
    size_t count;     // current number of lines
} line_batch_t;

typedef struct {
    float vertices[MAX_SPHERE_VERTICES];
    size_t vertex_count;
    size_t data_size;
} sphere_mesh_t;

// text rendering
typedef struct {
    GLuint tex, shader, vao, vbo;
    float verts[MAX_FONT_CHARS * 24];
    int count;
} font_t;

// planet path tracking
typedef struct {
    vec3 positions[MAX_SPACECRAFT * PATH_CAPACITY];
    int counts[MAX_SPACECRAFT];
    int capacity;
    int num_objects;
} craft_path_storage_t;
#endif


#endif
