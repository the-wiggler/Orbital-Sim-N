#ifndef CONFIG_H
#define CONFIG_H

#include <stdbool.h>

typedef struct {
    double time_step;
    int window_size_x, window_size_y;
    int screen_origin_x, screen_origin_y;
    double meters_per_pixel;
    int font_size;
    bool window_open;
    bool sim_running;
    double sim_time;
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
    double force_x;
    double force_y;
    double kinetic_energy;
} body_properties_t;

#endif