#include "spacecraft.h"
#include "../globals.h"
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include "../math/matrix.h"

void displayError(const char* title, const char* message);

// check and activate burns
void craft_checkBurnSchedule(const spacecraft_properties_t* sc, const int i, const body_properties_t* gb, const double sim_time) {
    // loop through all burns and check if any should be active
    bool burn_active = false;
    for (int j = 0; j < sc->num_burns[i]; j++) {
        // check if within the burn window
        if (sim_time >= sc->burn_properties[i][j].burn_start_time && sim_time < sc->burn_properties[i][j].burn_end_time && sc->fuel_mass[i] > 0) {
            sc->engine_on[i] = true;
            sc->throttle[i] = sc->burn_properties[i][j].throttle;

            // calculate attitude based on burn type
            quaternion_t final_attitude = {};
            const burn_properties_t* burn = &sc->burn_properties[i][j];
            const int target_id = burn->burn_target_id;

            if (burn->relative_burn_target.absolute) {
                // absolute: heading is in absolute space coordinates
                vec3 absolute_axis = {0.0, 0.0, 1.0};
                final_attitude = quaternionFromAxisAngle(absolute_axis, burn->burn_heading);
            } else if (burn->relative_burn_target.tangent) {
                // tangent: heading is relative to the velocity vector (tangent to orbit)
                const double rel_vel_x = sc->vel_x[i] - gb->vel_x[target_id];
                const double rel_vel_y = sc->vel_y[i] - gb->vel_y[target_id];
                const double rel_vel_z = sc->vel_z[i] - gb->vel_z[target_id];
                // first we create 3d velocity direction vector
                vec3 velocity_direction = {rel_vel_x, rel_vel_y, rel_vel_z};
                // default forward direction the spacecrafts local positive x axis, meaning the front of the craft
                vec3 default_forward = {1.0, 0.0, 0.0};
                // create quaternion to rotate from default forward to velocity direction
                quaternion_t base_rotation = quaternionFromTwoVectors(default_forward, velocity_direction);
                // apply burn heading offset
                if (burn->burn_heading != 0.0) {
                    // normalize velocity to get rotation axis
                    double vel_mag = sqrt(rel_vel_x * rel_vel_x + rel_vel_y * rel_vel_y + rel_vel_z * rel_vel_z);
                    vec3 rotation_axis = {rel_vel_x / vel_mag, rel_vel_y / vel_mag, rel_vel_z / vel_mag};
                    quaternion_t offset_rotation = quaternionFromAxisAngle(rotation_axis, burn->burn_heading);
                    // quaternion multiplication is right-to-left: offset * base means "first base, then offset"
                    final_attitude = quaternionMul(offset_rotation, base_rotation);
                } else {
                    final_attitude = base_rotation;
                }
            } else if (burn->relative_burn_target.normal) {
                displayError("ERROR", "Normal burns are currently unsupported. Please reconfigure. Sorry! :(");
            }
            else {
                displayError("ERROR", "Failed at determining burn type. If you see this you suck at coding lol");
            }

            sc->attitude[i] = final_attitude;
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
void craft_calculateGravForce(sim_properties_t* sim, const int craft_idx, const int body_idx) {
    const spacecraft_properties_t* sc = &sim->gs;
    const body_properties_t* bodies = &sim->gb;

    // calculate the distance between the spacecraft and the body
    const double delta_pos_x = bodies->pos_x[body_idx] - sc->pos_x[craft_idx];
    const double delta_pos_y = bodies->pos_y[body_idx] - sc->pos_y[craft_idx];
    const double delta_pos_z = bodies->pos_z[body_idx] - sc->pos_z[craft_idx];
    const double r_squared = delta_pos_x * delta_pos_x + delta_pos_y * delta_pos_y + delta_pos_z * delta_pos_z;

    // planet collision logic -- checks if craft is too close to the planet
    const double radius_squared = bodies->radius[body_idx] * bodies->radius[body_idx];
    if (r_squared < radius_squared) {
        sim->wp.sim_running = false;
        sim->wp.reset_sim = true;
        char err_txt[128];
        snprintf(err_txt, sizeof(err_txt), "Warning: %s has collided with %s\n\nResetting Simulation...", sc->names[craft_idx], bodies->names[body_idx]);
        displayError("PLANET COLLISION", err_txt);
        return;
    }

    // calculate the ship mass with the current amount of fuel in it
    sc->current_total_mass[craft_idx] = sc->fuel_mass[craft_idx] + sc->dry_mass[craft_idx];

    // force = (G * m1 * m2) * delta / r^3
    const double r = sqrt(r_squared);
    const double r_cubed = r_squared * r;
    const double force_factor = (G * sc->current_total_mass[craft_idx] * bodies->mass[body_idx]) / r_cubed;

    // apply the force to the craft
    sc->grav_force_x[craft_idx] += force_factor * delta_pos_x;
    sc->grav_force_y[craft_idx] += force_factor * delta_pos_y;
    sc->grav_force_z[craft_idx] += force_factor * delta_pos_z;
}

// updates the force
void craft_applyThrust(const spacecraft_properties_t* sc, const int i) {
    if (sc->engine_on[i] && sc->fuel_mass[i] > 0) {
        const double current_thrust = sc->thrust[i] * sc->throttle[i];

        // this defines the default vector position in which the engine applies the thrust relative to the spacecraft
        // in this case, setting x to 1 means that we are defining the "front" as the positive X direction relative to the craft.
        vec3 engine_thrust_direction = {1.0, 0.0, 0.0};

        // then we rotate the thrust direction by the spacecrafts defined attitude in space
        vec3 world_thrust = quaternionRotate(sc->attitude[i], engine_thrust_direction);

        // now we apply that thrust in the rotated direction in world space
        sc->grav_force_x[i] += current_thrust * world_thrust.x;
        sc->grav_force_y[i] += current_thrust * world_thrust.y;
        sc->grav_force_z[i] += current_thrust * world_thrust.z;
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
    sc->acc_z[i] = sc->grav_force_z[i] / sc->current_total_mass[i];

    // update position using current velocity and acceleration
    sc->pos_x[i] += (sc->vel_x[i] * dt) + (0.5 * sc->acc_x[i] * dt * dt);
    sc->pos_y[i] += (sc->vel_y[i] * dt) + (0.5 * sc->acc_y[i] * dt * dt);
    sc->pos_z[i] += (sc->vel_z[i] * dt) + (0.5 * sc->acc_z[i] * dt * dt);

    // update velocity using average of current and previous acceleration
    sc->vel_x[i] += 0.5 * (sc->acc_x[i] + sc->acc_x_prev[i]) * dt;
    sc->vel_y[i] += 0.5 * (sc->acc_y[i] + sc->acc_y_prev[i]) * dt;
    sc->vel_z[i] += 0.5 * (sc->acc_z[i] + sc->acc_z_prev[i]) * dt;
    sc->vel[i] = sqrt(sc->vel_x[i] * sc->vel_x[i] + sc->vel_y[i] * sc->vel_y[i] + sc->vel_z[i] * sc->vel_z[i]);

    // store current acceleration for next iteration
    sc->acc_x_prev[i] = sc->acc_x[i];
    sc->acc_y_prev[i] = sc->acc_y[i];
    sc->acc_z_prev[i] = sc->acc_z[i];
}

// adds a spacecraft to the spacecraft array
void craft_addSpacecraft(spacecraft_properties_t* sc, const char* name,
                        const double x_pos, const double y_pos, const double z_pos, const double x_vel, const double y_vel, const double z_vel,
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
    double* temp_pos_z = REALLOC_ARRAY(pos_z, double);
    quaternion_t* temp_attitude = REALLOC_ARRAY(attitude, quaternion_t);
    double* temp_vel_x = REALLOC_ARRAY(vel_x, double);
    double* temp_vel_y = REALLOC_ARRAY(vel_y, double);
    double* temp_vel_z = REALLOC_ARRAY(vel_z, double);
    double* temp_vel = REALLOC_ARRAY(vel, double);
    double* temp_rot_v = REALLOC_ARRAY(rotational_v, double);
    double* temp_momentum = REALLOC_ARRAY(momentum, double);
    double* temp_acc_x = REALLOC_ARRAY(acc_x, double);
    double* temp_acc_y = REALLOC_ARRAY(acc_y, double);
    double* temp_acc_z = REALLOC_ARRAY(acc_z, double);
    double* temp_acc_x_prev = REALLOC_ARRAY(acc_x_prev, double);
    double* temp_acc_y_prev = REALLOC_ARRAY(acc_y_prev, double);
    double* temp_acc_z_prev = REALLOC_ARRAY(acc_z_prev, double);
    double* temp_rot_a = REALLOC_ARRAY(rotational_a, double);
    double* temp_moi = REALLOC_ARRAY(moment_of_inertia, double);
    double* temp_grav_fx = REALLOC_ARRAY(grav_force_x, double);
    double* temp_grav_fy = REALLOC_ARRAY(grav_force_y, double);
    double* temp_grav_fz = REALLOC_ARRAY(grav_force_z, double);
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
    sc->pos_z = temp_pos_z;
    sc->attitude = temp_attitude;
    sc->vel_x = temp_vel_x;
    sc->vel_y = temp_vel_y;
    sc->vel_z = temp_vel_z;
    sc->vel = temp_vel;
    sc->rotational_v = temp_rot_v;
    sc->momentum = temp_momentum;
    sc->acc_x = temp_acc_x;
    sc->acc_y = temp_acc_y;
    sc->acc_z = temp_acc_z;
    sc->acc_x_prev = temp_acc_x_prev;
    sc->acc_y_prev = temp_acc_y_prev;
    sc->acc_z_prev = temp_acc_z_prev;
    sc->rotational_a = temp_rot_a;
    sc->moment_of_inertia = temp_moi;
    sc->grav_force_x = temp_grav_fx;
    sc->grav_force_y = temp_grav_fy;
    sc->grav_force_z = temp_grav_fz;
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
#ifdef WIN32
    strcpy_s(sc->names[idx], strlen(name) + 1, name);
#else
    strcpy(sc->names[idx], name);
#endif

    // initialize the new craft
    sc->pos_x[idx] = x_pos;
    sc->pos_y[idx] = y_pos;
    sc->pos_z[idx] = z_pos;
    sc->vel_x[idx] = x_vel;
    sc->vel_y[idx] = y_vel;
    sc->vel_z[idx] = z_vel;
    sc->vel[idx] = sqrt(x_vel * x_vel + y_vel * y_vel + z_vel * z_vel);
    sc->acc_x[idx] = 0.0;
    sc->acc_y[idx] = 0.0;
    sc->acc_z[idx] = 0.0;
    sc->acc_x_prev[idx] = 0.0;
    sc->acc_y_prev[idx] = 0.0;
    sc->acc_z_prev[idx] = 0.0;
    sc->grav_force_x[idx] = 0.0;
    sc->grav_force_y[idx] = 0.0;
    sc->grav_force_z[idx] = 0.0;
    vec3 start_axis = {0.0, 0.0, 1.0};
    sc->attitude[idx] = quaternionFromAxisAngle(start_axis, attitude);
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