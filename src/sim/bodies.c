#include "bodies.h"
#include "../globals.h"
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

void displayError(const char* title, const char* message);

// at the end of each sim loop, this function should be run to calculate the changes in
// the force values based on other parameters. for example, using F to find a based on m.
// i is the body that has the force applied to it, whilst j is the body applying force to i
void body_calculateGravForce(sim_properties_t* sim, const int i, const int j) {
    const body_properties_t* bodies = &sim->gb;

    // calculate the distance between the two bodies
    const double delta_pos_x = bodies->pos_x[j] - bodies->pos_x[i];
    const double delta_pos_y = bodies->pos_y[j] - bodies->pos_y[i];
    const double r_squared = delta_pos_x * delta_pos_x + delta_pos_y * delta_pos_y;

    // planet collision logic -- checks if planets are too close
    const double radius_squared = bodies->radius[i] * bodies->radius[i];
    if (r_squared < radius_squared) {
        sim->wp.sim_running = false;
        sim->wp.reset_sim = true;
        char err_txt[128];
        snprintf(err_txt, sizeof(err_txt), "Warning: %s has collided with %s\n\nResetting Simulation...", bodies->names[i], bodies->names[j]);
        displayError("PLANET COLLISION", err_txt);
        return;
    }

    // optimized version of the gravitation equation because computers dislike division
    // force = (G * m1 * m2) * delta / r^3
    const double r = sqrt(r_squared);
    const double r_cubed = r_squared * r;
    const double force_factor = (G * bodies->mass[i] * bodies->mass[j]) / r_cubed;

    const double force_x = force_factor * delta_pos_x;
    const double force_y = force_factor * delta_pos_y;

    // applies force to both bodies (since they apply the same force on one another)
    bodies->force_x[i] += force_x;
    bodies->force_y[i] += force_y;
    bodies->force_x[j] -= force_x;
    bodies->force_y[j] -= force_y;
}

// this calculates the changes of velocity and position based on the force values
// this uses a method called velocity verlet integration
void body_updateMotion(const body_properties_t* bodies, const int i, const double dt) {
    // calculate the current acceleration from the force on the object
    bodies->acc_x[i] = bodies->force_x[i] / bodies->mass[i];
    bodies->acc_y[i] = bodies->force_y[i] / bodies->mass[i];

    // update position using current velocity and acceleration
    bodies->pos_x[i] += (bodies->vel_x[i] * dt) + (0.5 * bodies->acc_x[i] * dt * dt);
    bodies->pos_y[i] += (bodies->vel_y[i] * dt) + (0.5 * bodies->acc_y[i] * dt * dt);

    // update velocity using average of current and previous acceleration
    // on first step, acc_prev will be zero, so we use current acceleration
    bodies->vel_x[i] += 0.5 * (bodies->acc_x[i] + bodies->acc_x_prev[i]) * dt;
    bodies->vel_y[i] += 0.5 * (bodies->acc_y[i] + bodies->acc_y_prev[i]) * dt;
    bodies->vel[i] = sqrt(bodies->vel_x[i] * bodies->vel_x[i] + bodies->vel_y[i] * bodies->vel_y[i]);

    // store current acceleration for next iteration
    bodies->acc_x_prev[i] = bodies->acc_x[i];
    bodies->acc_y_prev[i] = bodies->acc_y[i];
}

// calculates the kinetic energy of a target body
void body_calculateKineticEnergy(const body_properties_t* bodies, const int i) {
    // calculate kinetic energy (0.5mv^2)
    bodies->kinetic_energy[i] = 0.5 * bodies->mass[i] * bodies->vel[i] * bodies->vel[i];
}

// calculates the size (in pixels) that the planet should appear on the screen based on its mass
float body_calculateVisualRadius(const body_properties_t* bodies, const int i, const window_params_t wp) {
    float r = (float)bodies->radius[i] / (float)wp.meters_per_pixel;
    bodies->pixel_radius[i] = r;
    return r;
}

// function to add a new body to the system
void body_addOrbitalBody(body_properties_t* gb, const char* name, const double mass, const double radius, const double x_pos, const double y_pos, const double x_vel, const double y_vel) {
    int new_size = gb->count + 1;

    // reallocate all arrays
    char** temp_names = (char**)realloc(gb->names, new_size * sizeof(char*));
    double* temp_mass = (double*)realloc(gb->mass, new_size * sizeof(double));
    double* temp_radius = (double*)realloc(gb->radius, new_size * sizeof(double));
    float* temp_pixel_radius = (float*)realloc(gb->pixel_radius, new_size * sizeof(float));
    double* temp_pos_x = (double*)realloc(gb->pos_x, new_size * sizeof(double));
    double* temp_pos_y = (double*)realloc(gb->pos_y, new_size * sizeof(double));
    double* temp_vel_x = (double*)realloc(gb->vel_x, new_size * sizeof(double));
    double* temp_vel_y = (double*)realloc(gb->vel_y, new_size * sizeof(double));
    double* temp_vel = (double*)realloc(gb->vel, new_size * sizeof(double));
    double* temp_acc_x = (double*)realloc(gb->acc_x, new_size * sizeof(double));
    double* temp_acc_y = (double*)realloc(gb->acc_y, new_size * sizeof(double));
    double* temp_acc_x_prev = (double*)realloc(gb->acc_x_prev, new_size * sizeof(double));
    double* temp_acc_y_prev = (double*)realloc(gb->acc_y_prev, new_size * sizeof(double));
    double* temp_force_x = (double*)realloc(gb->force_x, new_size * sizeof(double));
    double* temp_force_y = (double*)realloc(gb->force_y, new_size * sizeof(double));
    double* temp_kinetic = (double*)realloc(gb->kinetic_energy, new_size * sizeof(double));
    float** temp_cached_body_coords_x = (float**)realloc(gb->cached_body_coords_x, new_size * sizeof(float*));
    float** temp_cached_body_coords_y = (float**)realloc(gb->cached_body_coords_y, new_size * sizeof(float*));

    if (!temp_names || !temp_mass || !temp_radius || !temp_pixel_radius ||
        !temp_pos_x || !temp_pos_y ||
        !temp_vel_x || !temp_vel_y || !temp_vel || !temp_acc_x || !temp_acc_y ||
        !temp_acc_x_prev || !temp_acc_y_prev || !temp_force_x || !temp_force_y || !temp_kinetic ||
        !temp_cached_body_coords_x || !temp_cached_body_coords_y) {
        displayError("ERROR", "Error: Failed to allocate memory for new body\n");
        return;
    }

    // allocate cache arrays for the new body
    temp_cached_body_coords_x[gb->count] = (float*)calloc(PATH_CACHE_LENGTH, sizeof(float));
    temp_cached_body_coords_y[gb->count] = (float*)calloc(PATH_CACHE_LENGTH, sizeof(float));
    if (!temp_cached_body_coords_x[gb->count] || !temp_cached_body_coords_y[gb->count]) {
        displayError("ERROR", "Error: Failed to allocate memory for body cache arrays\n");
        return;
    }

    gb->names = temp_names;
    gb->mass = temp_mass;
    gb->radius = temp_radius;
    gb->pixel_radius = temp_pixel_radius;
    gb->pos_x = temp_pos_x;
    gb->pos_y = temp_pos_y;
    gb->vel_x = temp_vel_x;
    gb->vel_y = temp_vel_y;
    gb->vel = temp_vel;
    gb->acc_x = temp_acc_x;
    gb->acc_y = temp_acc_y;
    gb->acc_x_prev = temp_acc_x_prev;
    gb->acc_y_prev = temp_acc_y_prev;
    gb->force_x = temp_force_x;
    gb->force_y = temp_force_y;
    gb->kinetic_energy = temp_kinetic;
    gb->cached_body_coords_x = temp_cached_body_coords_x;
    gb->cached_body_coords_y = temp_cached_body_coords_y;

    int idx = gb->count;

    // allocate memory for the name and copy it
    gb->names[idx] = (char*)malloc(strlen(name) + 1);
    if (gb->names[idx] == NULL) {
        displayError("ERROR", "Error: Failed to allocate memory for body name\n");
        return;
    }
#ifdef WIN32
    strcpy_s(gb->names[idx], strlen(name) + 1, name);
#else
    strcpy(gb->names[idx], name);
#endif

    // initialize the new body
    gb->mass[idx] = mass;
    gb->radius[idx] = radius;
    gb->pixel_radius[idx] = 0.0f;
    gb->pos_x[idx] = x_pos;
    gb->pos_y[idx] = y_pos;
    gb->vel_x[idx] = x_vel;
    gb->vel_y[idx] = y_vel;
    gb->vel[idx] = sqrt(x_vel * x_vel + y_vel * y_vel);
    gb->acc_x[idx] = 0.0;
    gb->acc_y[idx] = 0.0;
    gb->acc_x_prev[idx] = 0.0;
    gb->acc_y_prev[idx] = 0.0;
    gb->force_x[idx] = 0.0;
    gb->force_y[idx] = 0.0;
    gb->kinetic_energy[idx] = 0.5 * mass * gb->vel[idx] * gb->vel[idx];

    // increment the body count
    gb->count++;
}
