#include "spacecraft.h"
#include "../globals.h"
#include <math.h>
#include <string.h>
#include <stdio.h>

void displayError(const char* title, const char* message);

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
                        const burn_properties_t* burns, const int num_burns) {

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

    // initialize burn schedule with provided burns
    sc->num_burns[idx] = num_burns;
    if (num_burns > 0) {
        sc->burn_properties[idx] = (burn_properties_t*)malloc(num_burns * sizeof(burn_properties_t));
        for (int i = 0; i < num_burns; i++) {
            sc->burn_properties[idx][i] = burns[i];
        }
    } else {
        sc->burn_properties[idx] = NULL;
    }

    // increment the craft count
    sc->count++;
}

// transform the craft coordinates in meters to pixel coordinates on the screen
void craft_transformCoordinates(const spacecraft_properties_t* sc, const int i, const window_params_t wp) {
    sc->pixel_coordinates_x[i] = wp.screen_origin_x + (float)(sc->pos_x[i] / wp.meters_per_pixel);
    sc->pixel_coordinates_y[i] = wp.screen_origin_y - (float)(sc->pos_y[i] / wp.meters_per_pixel); // this is negative because the SDL origin is in the top left, so positive y is 'down'
}
