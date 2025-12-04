////////////////////////////////////////////////////////////////////////////////////////////////////
// SIMULATION CALCULATION FUNCTIONS
////////////////////////////////////////////////////////////////////////////////////////////////////
#include "sim_calculations.h"
#include "config.h"
#include "sdl_elements.h"
#include <math.h>
#include <string.h>
#include <immintrin.h>
#include <stdio.h>

const double G = 6.67430E-11;

////////////////////////////////////////////////////////////////////////////////////////////////////
// ORBITAL BODY CALCULATIONS
////////////////////////////////////////////////////////////////////////////////////////////////////
// at the end of each sim loop, this function should be run to calculate the changes in
// the force values based on other parameters. for example, using F to find a based on m.
// i is the body that has the force applied to it, whilst j is the body applying force to i
void body_calculateGravForce(const body_properties_t* bodies, const int i, const int j, window_params_t* wp) {
    // calculate the distance between the two bodies
    const double delta_pos_x = bodies->pos_x[j] - bodies->pos_x[i];
    const double delta_pos_y = bodies->pos_y[j] - bodies->pos_y[i];
    const double r_squared = delta_pos_x * delta_pos_x + delta_pos_y * delta_pos_y;

    // planet collision logic -- checks if planets are too close
    const double radius_squared = bodies->radius[i] * bodies->radius[i];
    if (r_squared < radius_squared) {
        wp->sim_running = false;
        wp->reset_sim = true;
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

// calculate total system energy for all bodies
// this avoids double-counting by only calculating each pair interaction once
double calculateTotalSystemEnergy(const body_properties_t* gb, const spacecraft_properties_t* sc) {
    double total_kinetic = 0.0;
    double total_potential = 0.0;

    // calculate kinetic energy for all bodies
    for (int i = 0; i < gb->count; i++) {
        total_kinetic += 0.5 * gb->mass[i] * gb->vel[i] * gb->vel[i];
    }

    // calculate kinetic energy for all spacecraft
    for (int i = 0; i < sc->count; i++) {
        total_kinetic += 0.5 * sc->current_total_mass[i] * sc->vel[i] * sc->vel[i];
    }

    // calculate potential energy between all body pairs
    for (int i = 0; i < gb->count; i++) {
        for (int j = i + 1; j < gb->count; j++) {
            double delta_x = gb->pos_x[j] - gb->pos_x[i];
            double delta_y = gb->pos_y[j] - gb->pos_y[i];
            double r = sqrt(delta_x * delta_x + delta_y * delta_y);
            if (r > 0) {  // prevent divide by zero
                total_potential += -(G * gb->mass[i] * gb->mass[j]) / r;
            }
        }
    }

    // calculate potential energy between spacecraft and bodies
    for (int i = 0; i < sc->count; i++) {
        for (int j = 0; j < gb->count; j++) {
            double delta_x = gb->pos_x[j] - sc->pos_x[i];
            double delta_y = gb->pos_y[j] - sc->pos_y[i];
            double r = sqrt(delta_x * delta_x + delta_y * delta_y);
            if (r > 0) {  // prevent divide by zero
                total_potential += -(G * sc->current_total_mass[i] * gb->mass[j]) / r;
            }
        }
    }

    return total_kinetic + total_potential;
}

// transforms spacial coordinates (for example, in meters) to pixel coordinates
void body_transformCoordinates(const body_properties_t* bodies, const int i, const window_params_t wp) {
    bodies->pixel_coordinates_x[i] = wp.screen_origin_x + ((float)bodies->pos_x[i] / (float)wp.meters_per_pixel);
    bodies->pixel_coordinates_y[i] = wp.screen_origin_y - ((float)bodies->pos_y[i] / (float)wp.meters_per_pixel); // this is negative because the SDL origin is in the top left, so positive y is 'down'
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
    float* temp_pixel_x = (float*)realloc(gb->pixel_coordinates_x, new_size * sizeof(float));
    float* temp_pixel_y = (float*)realloc(gb->pixel_coordinates_y, new_size * sizeof(float));
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
        !temp_pos_x || !temp_pos_y || !temp_pixel_x || !temp_pixel_y ||
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
    gb->pixel_coordinates_x = temp_pixel_x;
    gb->pixel_coordinates_y = temp_pixel_y;
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
    strcpy_s(gb->names[idx], strlen(name) + 1, name);

    // initialize the new body
    gb->mass[idx] = mass;
    gb->radius[idx] = radius;
    gb->pixel_radius[idx] = 0.0f;
    gb->pos_x[idx] = x_pos;
    gb->pos_y[idx] = y_pos;
    gb->pixel_coordinates_x[idx] = 0.0f;
    gb->pixel_coordinates_y[idx] = 0.0f;
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

////////////////////////////////////////////////////////////////////////////////////////////////////
// SPACECRAFT CALCULATIONS
////////////////////////////////////////////////////////////////////////////////////////////////////
// calculates the force applied on a spacecraft by a specific body
// (this should be used before any other force functions)
void craft_calculateGravForce(const spacecraft_properties_t* sc, const int i, const body_properties_t* bodies, const int j, window_params_t* wp) {
    // calculate the distance between the spacecraft and the body
    const double delta_pos_x = bodies->pos_x[j] - sc->pos_x[i];
    const double delta_pos_y = bodies->pos_y[j] - sc->pos_y[i];
    const double r_squared = delta_pos_x * delta_pos_x + delta_pos_y * delta_pos_y;

    // planet collision logic -- checks if craft is too close to the planet
    const double radius_squared = bodies->radius[j] * bodies->radius[j];
    if (r_squared < radius_squared) {
        wp->sim_running = false;
        wp->reset_sim = true;
        char err_txt[128];
        snprintf(err_txt, sizeof(err_txt), "Warning: %s has collided with %s\n\nResetting Simulation...", sc->names[i], bodies->names[j]);
        displayError("PLANET COLLISION", err_txt);
        return;
    }

    // calculate the ship mass with the current amount of fuel in it
    sc->current_total_mass[i] = sc->fuel_mass[i] + sc->dry_mass[i];

    // prevent division by zero
    const double min_distance_squared = 1.0;
    if (r_squared < min_distance_squared) {
        return;  // skip calculation if spacecraft is too close to body
    }

    // force = (G * m1 * m2) * delta / r^3
    const double r = sqrt(r_squared);
    const double r_cubed = r_squared * r;
    const double force_factor = (G * sc->current_total_mass[i] * bodies->mass[j]) / r_cubed;

    sc->grav_force_x[i] += force_factor * delta_pos_x;
    sc->grav_force_y[i] += force_factor * delta_pos_y;
}

// updates the force
void craft_applyThrust(const spacecraft_properties_t* sc, const int i) {
    if (sc->engine_on[i] && sc->fuel_mass[i] > 0) {
        const double current_thrust = sc->thrust[i] * sc->throttle[i];
        sc->grav_force_x[i] += current_thrust * cos(sc->attitude[i]);
        sc->grav_force_y[i] += current_thrust * sin(sc->attitude[i]);
    }
}

// check and activate burns
void craft_checkBurnSchedule(const spacecraft_properties_t* sc, const int i, const double sim_time) {
    // loop through all burns and check if any should be active
    bool burn_active = false;
    for (int j = 0; j < sc->num_burns[i]; j++) {
        // check if within the burn window
        if (sim_time >= sc->burn_properties[i][j].burn_start_time &&
            sim_time < sc->burn_properties[i][j].burn_end_time &&
            sc->fuel_mass[i] > 0) {
            sc->engine_on[i] = true;
            sc->throttle[i] = sc->burn_properties[i][j].throttle;
            sc->attitude[i] = sc->burn_properties[i][j].burn_heading;
            burn_active = true;
            break; // only execute one burn at a time
        }
    }

    // if no burn is active, turn off the engine
    if (!burn_active) {
        sc->engine_on[i] = false;
        sc->throttle[i] = 0.0;
    }
}


void craft_consumeFuel(const spacecraft_properties_t* sc, const int i, const double dt) {
    if (sc->engine_on[i] && sc->fuel_mass[i] > 0) {
        double fuel_consumed = sc->mass_flow_rate[i] * sc->throttle[i] * dt;

        if (fuel_consumed > sc->fuel_mass[i]) {
            fuel_consumed = sc->fuel_mass[i];
            sc->engine_on[i] = false;
        }

        sc->fuel_mass[i] -= fuel_consumed;
        sc->current_total_mass[i] = sc->dry_mass[i] + sc->fuel_mass[i];
    }
}


// updates the motion of the spacecraft based on the force currently applied to it
// uses velocity verlet integration
void craft_updateMotion(const spacecraft_properties_t* sc, const int i, const double dt) {
    // calculate the current acceleration from the force on the object
    sc->acc_x[i] = sc->grav_force_x[i] / sc->current_total_mass[i];
    sc->acc_y[i] = sc->grav_force_y[i] / sc->current_total_mass[i];

    // update position using current velocity and acceleration
    sc->pos_x[i] += (sc->vel_x[i] * dt) + (0.5 * sc->acc_x[i] * dt * dt);
    sc->pos_y[i] += (sc->vel_y[i] * dt) + (0.5 * sc->acc_y[i] * dt * dt);

    // update velocity using average of current and previous acceleration
    sc->vel_x[i] += 0.5 * (sc->acc_x[i] + sc->acc_x_prev[i]) * dt;
    sc->vel_y[i] += 0.5 * (sc->acc_y[i] + sc->acc_y_prev[i]) * dt;
    sc->vel[i] = sqrt(sc->vel_x[i] * sc->vel_x[i] + sc->vel_y[i] * sc->vel_y[i]);

    // store current acceleration for next iteration
    sc->acc_x_prev[i] = sc->acc_x[i];
    sc->acc_y_prev[i] = sc->acc_y[i];
}

// adds a spacecraft to the spacecraft array
void craft_addSpacecraft(spacecraft_properties_t* sc, const char* name,
                        const double x_pos, const double y_pos, const double x_vel, const double y_vel,
                        const double dry_mass, const double fuel_mass, const double thrust,
                        const double specific_impulse, const double mass_flow_rate,
                        const double attitude, const double moment_of_inertia,
                        const double nozzle_gimbal_range,
                        const double burn_start_time, const double burn_duration,
                        const double burn_heading, const double burn_throttle) {

    int new_size = sc->count + 1;

    // reallocate all arrays
    #define REALLOC_ARRAY(arr, type) (type*)realloc(sc->arr, new_size * sizeof(type))

    char** temp_names = REALLOC_ARRAY(names, char*);
    double* temp_current_mass = REALLOC_ARRAY(current_total_mass, double);
    double* temp_dry_mass = REALLOC_ARRAY(dry_mass, double);
    double* temp_fuel_mass = REALLOC_ARRAY(fuel_mass, double);
    double* temp_pos_x = REALLOC_ARRAY(pos_x, double);
    double* temp_pos_y = REALLOC_ARRAY(pos_y, double);
    float* temp_pixel_x = REALLOC_ARRAY(pixel_coordinates_x, float);
    float* temp_pixel_y = REALLOC_ARRAY(pixel_coordinates_y, float);
    double* temp_attitude = REALLOC_ARRAY(attitude, double);
    double* temp_vel_x = REALLOC_ARRAY(vel_x, double);
    double* temp_vel_y = REALLOC_ARRAY(vel_y, double);
    double* temp_vel = REALLOC_ARRAY(vel, double);
    double* temp_rot_v = REALLOC_ARRAY(rotational_v, double);
    double* temp_momentum = REALLOC_ARRAY(momentum, double);
    double* temp_acc_x = REALLOC_ARRAY(acc_x, double);
    double* temp_acc_y = REALLOC_ARRAY(acc_y, double);
    double* temp_acc_x_prev = REALLOC_ARRAY(acc_x_prev, double);
    double* temp_acc_y_prev = REALLOC_ARRAY(acc_y_prev, double);
    double* temp_rot_a = REALLOC_ARRAY(rotational_a, double);
    double* temp_moi = REALLOC_ARRAY(moment_of_inertia, double);
    double* temp_grav_fx = REALLOC_ARRAY(grav_force_x, double);
    double* temp_grav_fy = REALLOC_ARRAY(grav_force_y, double);
    double* temp_torque = REALLOC_ARRAY(torque, double);
    double* temp_thrust = REALLOC_ARRAY(thrust, double);
    double* temp_mfr = REALLOC_ARRAY(mass_flow_rate, double);
    double* temp_isp = REALLOC_ARRAY(specific_impulse, double);
    double* temp_throttle = REALLOC_ARRAY(throttle, double);
    double* temp_gimbal = REALLOC_ARRAY(nozzle_gimbal_range, double);
    double* temp_nozzle_v = REALLOC_ARRAY(nozzle_velocity, double);
    bool* temp_engine = REALLOC_ARRAY(engine_on, bool);
    int* temp_num_burns = REALLOC_ARRAY(num_burns, int);
    burn_properties_t** temp_burns = REALLOC_ARRAY(burn_properties, burn_properties_t*);

    #undef REALLOC_ARRAY

    sc->names = temp_names;
    sc->current_total_mass = temp_current_mass;
    sc->dry_mass = temp_dry_mass;
    sc->fuel_mass = temp_fuel_mass;
    sc->pos_x = temp_pos_x;
    sc->pos_y = temp_pos_y;
    sc->pixel_coordinates_x = temp_pixel_x;
    sc->pixel_coordinates_y = temp_pixel_y;
    sc->attitude = temp_attitude;
    sc->vel_x = temp_vel_x;
    sc->vel_y = temp_vel_y;
    sc->vel = temp_vel;
    sc->rotational_v = temp_rot_v;
    sc->momentum = temp_momentum;
    sc->acc_x = temp_acc_x;
    sc->acc_y = temp_acc_y;
    sc->acc_x_prev = temp_acc_x_prev;
    sc->acc_y_prev = temp_acc_y_prev;
    sc->rotational_a = temp_rot_a;
    sc->moment_of_inertia = temp_moi;
    sc->grav_force_x = temp_grav_fx;
    sc->grav_force_y = temp_grav_fy;
    sc->torque = temp_torque;
    sc->thrust = temp_thrust;
    sc->mass_flow_rate = temp_mfr;
    sc->specific_impulse = temp_isp;
    sc->throttle = temp_throttle;
    sc->nozzle_gimbal_range = temp_gimbal;
    sc->nozzle_velocity = temp_nozzle_v;
    sc->engine_on = temp_engine;
    sc->num_burns = temp_num_burns;
    sc->burn_properties = temp_burns;

    int idx = sc->count;

    // allocate memory for the name and copy it
    sc->names[idx] = (char*)malloc(strlen(name) + 1);
    strcpy_s(sc->names[idx], strlen(name) + 1, name);

    // initialize the new craft
    sc->pos_x[idx] = x_pos;
    sc->pos_y[idx] = y_pos;
    sc->vel_x[idx] = x_vel;
    sc->vel_y[idx] = y_vel;
    sc->vel[idx] = sqrt(x_vel * x_vel + y_vel * y_vel);
    sc->acc_x[idx] = 0.0;
    sc->acc_y[idx] = 0.0;
    sc->acc_x_prev[idx] = 0.0;
    sc->acc_y_prev[idx] = 0.0;
    sc->grav_force_x[idx] = 0.0;
    sc->grav_force_y[idx] = 0.0;
    sc->pixel_coordinates_x[idx] = 0.0f;
    sc->pixel_coordinates_y[idx] = 0.0f;
    sc->attitude[idx] = attitude;
    sc->dry_mass[idx] = dry_mass;
    sc->fuel_mass[idx] = fuel_mass;
    sc->current_total_mass[idx] = dry_mass + fuel_mass;
    sc->mass_flow_rate[idx] = mass_flow_rate;
    sc->thrust[idx] = thrust;
    sc->specific_impulse[idx] = specific_impulse;
    sc->throttle[idx] = 0.0;
    sc->engine_on[idx] = false;
    sc->rotational_v[idx] = 0.0;
    sc->momentum[idx] = 0.0;
    sc->rotational_a[idx] = 0.0;
    sc->moment_of_inertia[idx] = moment_of_inertia;
    sc->torque[idx] = 0.0;
    sc->nozzle_gimbal_range[idx] = nozzle_gimbal_range;
    sc->nozzle_velocity[idx] = 0.0;

    // initialize burn schedule with a single burn
    sc->num_burns[idx] = 1;
    sc->burn_properties[idx] = (burn_properties_t*)malloc(sizeof(burn_properties_t));
    sc->burn_properties[idx][0].burn_start_time = burn_start_time;
    sc->burn_properties[idx][0].burn_end_time = burn_start_time + burn_duration;
    sc->burn_properties[idx][0].throttle = burn_throttle;
    sc->burn_properties[idx][0].burn_heading = burn_heading;

    // increment the craft count
    sc->count++;
}

// transform the craft coordinates in meters to pixel coordinates on the screen
void craft_transformCoordinates(const spacecraft_properties_t* sc, const int i, const window_params_t wp) {
    sc->pixel_coordinates_x[i] = wp.screen_origin_x + (float)(sc->pos_x[i] / wp.meters_per_pixel);
    sc->pixel_coordinates_y[i] = wp.screen_origin_y - (float)(sc->pos_y[i] / wp.meters_per_pixel); // this is negative because the SDL origin is in the top left, so positive y is 'down'
}


////////////////////////////////////////////////////////////////////////////////////////////////////
// LOGIC FUNCTIONS
////////////////////////////////////////////////////////////////////////////////////////////////////
// reset the simulation by removing all bodies from the system
void resetSim(window_params_t* wp, body_properties_t* gb, spacecraft_properties_t* sc) {
    // reset simulation time to 0
    wp->sim_time = 0;
    // change sim reset flag back to false
    wp->reset_sim = false;

    // free all bodies from memory
    if (gb != NULL) {
        for (int i = 0; i < gb->count; i++) {
            free(gb->names[i]);
        }

        // free cache arrays for each body
        if (gb->cached_body_coords_x != NULL) {
            free(gb->cached_body_coords_x);
        }
        if (gb->cached_body_coords_y != NULL) {
            free(gb->cached_body_coords_y);
        }

        free(gb->names);
        free(gb->mass);
        free(gb->radius);
        free(gb->pixel_radius);
        free(gb->pos_x);
        free(gb->pos_y);
        free(gb->pixel_coordinates_x);
        free(gb->pixel_coordinates_y);
        free(gb->vel_x);
        free(gb->vel_y);
        free(gb->vel);
        free(gb->acc_x);
        free(gb->acc_y);
        free(gb->acc_x_prev);
        free(gb->acc_y_prev);
        free(gb->force_x);
        free(gb->force_y);
        free(gb->kinetic_energy);

        gb->count = 0;
        gb->names = NULL;
        gb->mass = NULL;
        gb->radius = NULL;
        gb->pixel_radius = NULL;
        gb->pos_x = NULL;
        gb->pos_y = NULL;
        gb->pixel_coordinates_x = NULL;
        gb->pixel_coordinates_y = NULL;
        gb->vel_x = NULL;
        gb->vel_y = NULL;
        gb->vel = NULL;
        gb->acc_x = NULL;
        gb->acc_y = NULL;
        gb->acc_x_prev = NULL;
        gb->acc_y_prev = NULL;
        gb->force_x = NULL;
        gb->force_y = NULL;
        gb->kinetic_energy = NULL;
        gb->cached_body_coords_x = NULL;
        gb->cached_body_coords_y = NULL;
        gb->cache_counter = 0;
        gb->cache_valid_count = 0;
    }

    // free all craft from memory
    if (sc != NULL) {
        for (int i = 0; i < sc->count; i++) {
            free(sc->names[i]);
            free(sc->burn_properties[i]);
        }
        free(sc->names);
        free(sc->current_total_mass);
        free(sc->dry_mass);
        free(sc->fuel_mass);
        free(sc->pos_x);
        free(sc->pos_y);
        free(sc->pixel_coordinates_x);
        free(sc->pixel_coordinates_y);
        free(sc->attitude);
        free(sc->vel_x);
        free(sc->vel_y);
        free(sc->vel);
        free(sc->rotational_v);
        free(sc->momentum);
        free(sc->acc_x);
        free(sc->acc_y);
        free(sc->acc_x_prev);
        free(sc->acc_y_prev);
        free(sc->rotational_a);
        free(sc->moment_of_inertia);
        free(sc->grav_force_x);
        free(sc->grav_force_y);
        free(sc->torque);
        free(sc->thrust);
        free(sc->mass_flow_rate);
        free(sc->specific_impulse);
        free(sc->throttle);
        free(sc->nozzle_gimbal_range);
        free(sc->nozzle_velocity);
        free(sc->engine_on);
        free(sc->num_burns);
        free(sc->burn_properties);

        sc->count = 0;
        sc->names = NULL;
        sc->current_total_mass = NULL;
        sc->dry_mass = NULL;
        sc->fuel_mass = NULL;
        sc->pos_x = NULL;
        sc->pos_y = NULL;
        sc->pixel_coordinates_x = NULL;
        sc->pixel_coordinates_y = NULL;
        sc->attitude = NULL;
        sc->vel_x = NULL;
        sc->vel_y = NULL;
        sc->vel = NULL;
        sc->rotational_v = NULL;
        sc->momentum = NULL;
        sc->acc_x = NULL;
        sc->acc_y = NULL;
        sc->acc_x_prev = NULL;
        sc->acc_y_prev = NULL;
        sc->rotational_a = NULL;
        sc->moment_of_inertia = NULL;
        sc->grav_force_x = NULL;
        sc->grav_force_y = NULL;
        sc->torque = NULL;
        sc->thrust = NULL;
        sc->mass_flow_rate = NULL;
        sc->specific_impulse = NULL;
        sc->throttle = NULL;
        sc->nozzle_gimbal_range = NULL;
        sc->nozzle_velocity = NULL;
        sc->engine_on = NULL;
        sc->num_burns = NULL;
        sc->burn_properties = NULL;
    }
}

// calculate the optimum velocity for an object to orbit a given body based on the orbit radius
// (function does not exist yet)

///////////////////////////////////////////////////////////////////////////////////////////////////////
// MAIN CALCULATION LOOP - THIS LOOP IS PLACED IN THE MAIN LOOP, AND DOES ALL NECESSARY CALCULATIONS //
///////////////////////////////////////////////////////////////////////////////////////////////////////
void runCalculations(const body_properties_t* gb, const spacecraft_properties_t* sc, window_params_t* wp) {
    if (wp->sim_running) {
        ////////////////////////////////////////////////////////////////
        // calculate forces between all body pairs
        ////////////////////////////////////////////////////////////////
        if (gb != NULL && gb->count > 0) {
            // initialize forces to zero to re-calculate them
            for (int i = 0; i < gb->count; i++) {
                gb->force_x[i] = 0;
                gb->force_y[i] = 0;
            }

            // loop through every body and add the resultant force to the subject body force vector
            for (int i = 0; i < gb->count; i++) {
                for (int j = i + 1; j < gb->count; j++) {
                    body_calculateGravForce(gb, i, j, wp);
                }
            }

            // calculate kinetic energy for each body
            for (int i = 0; i < gb->count; i++) {
                body_calculateKineticEnergy(gb, i);
            }

            // update the motion for each body and draw
            for (int i = 0; i < gb->count; i++) {
                // updates the kinematic properties of each body (velocity, acceleration, position, etc.)
                body_updateMotion(gb, i, wp->time_step);
                // transform real-space coordinate to pixel coordinates on screen (scaling)
                body_transformCoordinates(gb, i, *wp);
            }
        }

        ////////////////////////////////////////////////////////////////
        // calculate forces between spacecraft and bodies
        ////////////////////////////////////////////////////////////////
        if (sc != NULL && sc->count > 0 && gb != NULL && gb->count > 0) {
            for (int i = 0; i < sc->count; i++) {
                sc->grav_force_x[i] = 0;
                sc->grav_force_y[i] = 0;

                // check if burn should be active based on simulation time
                craft_checkBurnSchedule(sc, i, wp->sim_time);

                // loop through all bodies and calculate gravitational forces on spacecraft
                for (int j = 0; j < gb->count; j++) {
                    craft_calculateGravForce(sc, i, gb, j, wp);
                }

                // apply thrust first, then consume fuel
                craft_applyThrust(sc, i);
                craft_consumeFuel(sc, i, wp->time_step);
            }

            // update the motion for each craft and draw
            for (int i = 0; i < sc->count; i++) {
                // updates the kinematic properties of each body (velocity, acceleration, position, etc.)
                craft_updateMotion(sc, i, wp->time_step);
                // transform real-space coordinate to pixel coordinates on screen (scaling)
                craft_transformCoordinates(sc, i, *wp);
            }
        }

        // increment the time based on the time step
        if (gb != NULL && gb->count > 0) {
            wp->sim_time += wp->time_step;
        }
    }
}