#include "spacecraft.h"
#include "../globals.h"
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include "../math/matrix.h"

void displayError(const char* title, const char* message);

// check and activate burns
void craft_checkBurnSchedule(spacecraft_t* craft, const body_properties_t* gb, const double sim_time) {
    bool burn_active = false;
    for (int j = 0; j < craft->num_burns; j++) {
        burn_properties_t* burn = &craft->burn_properties[j];
        // check if within the burn window
        if (sim_time >= burn->burn_start_time && sim_time < burn->burn_end_time && craft->fuel_mass > 0) {
            craft->engine_on = true;
            craft->throttle = burn->throttle;

            // calculate attitude based on burn type
            quaternion_t final_attitude = {0};
            int target_id = burn->burn_target_id;
            body_t* target = &gb->bodies[target_id];

            if (burn->relative_burn_target.absolute) {
                // absolute: heading is in absolute space coordinates
                vec3 absolute_axis = {0.0, 0.0, 1.0};
                final_attitude = quaternionFromAxisAngle(absolute_axis, burn->burn_heading);
            } else if (burn->relative_burn_target.tangent) {
                // tangent: heading is relative to the velocity vector
                vec3 rel_vel = vec3_sub(craft->vel, target->vel);
                vec3 default_forward = {0.0, 1.0, 0.0};
                quaternion_t base_rotation = quaternionFromTwoVectors(default_forward, rel_vel);

                if (burn->burn_heading != 0.0) {
                    vec3 rotation_axis = vec3_normalize(rel_vel);
                    quaternion_t offset_rotation = quaternionFromAxisAngle(rotation_axis, burn->burn_heading);
                    final_attitude = quaternionMul(offset_rotation, base_rotation);
                } else {
                    final_attitude = base_rotation;
                }
            } else if (burn->relative_burn_target.normal) {
                // normal: heading is perpendicular to the orbital plane
                vec3 rel_pos = vec3_sub(craft->pos, target->pos);
                vec3 rel_vel = vec3_sub(craft->vel, target->vel);
                vec3 normal_direction = cross_product_vec3(rel_pos, rel_vel);
                vec3 default_forward = {0.0, 1.0, 0.0};
                quaternion_t base_rotation = quaternionFromTwoVectors(default_forward, normal_direction);

                if (burn->burn_heading != 0.0) {
                    vec3 rotation_axis = vec3_normalize(normal_direction);
                    quaternion_t offset_rotation = quaternionFromAxisAngle(rotation_axis, burn->burn_heading);
                    final_attitude = quaternionMul(offset_rotation, base_rotation);
                } else {
                    final_attitude = base_rotation;
                }
            } else {
                displayError("ERROR", "Failed at determining burn type. If you see this you suck at coding lol");
            }

            craft->attitude = final_attitude;
            burn_active = true;
            break; // only execute one burn at a time
        }
    }

    // if no burn is active, turn off the engine
    if (!burn_active) {
        craft->engine_on = false;
        craft->throttle = 0.0;
    }
}

// calculates the force applied on a spacecraft by a specific body
void craft_calculateGravForce(sim_properties_t* sim, const int craft_idx, const int body_idx) {
    spacecraft_t* craft = &sim->gs.spacecraft[craft_idx];
    body_t* body = &sim->gb.bodies[body_idx];

    // calculate the distance between the spacecraft and the body
    vec3 delta_pos = vec3_sub(body->pos, craft->pos);
    double r_squared = vec3_mag_sq(delta_pos);
    double r = sqrt(r_squared);

    // planet collision logic
    if (r_squared < body->radius * body->radius) {
        sim->wp.sim_running = false;
        sim->wp.reset_sim = true;
        char err_txt[128];
        snprintf(err_txt, sizeof(err_txt), "Warning: %s has collided with %s\n\nResetting Simulation...", craft->name, body->name);
        displayError("PLANET COLLISION", err_txt);
        return;
    }

    // calculate the ship mass with the current amount of fuel
    craft->current_total_mass = craft->fuel_mass + craft->dry_mass;

    // force = (G * m1 * m2) * delta / r^3
    double r_cubed = r_squared * r;
    double force_factor = (G * craft->current_total_mass * body->mass) / r_cubed;

    // apply the force to the craft
    vec3 force = vec3_scale(delta_pos, force_factor);
    craft->grav_force = vec3_add(craft->grav_force, force);
}

// applies thrust force based on current attitude
void craft_applyThrust(spacecraft_t* craft) {
    if (craft->engine_on && craft->fuel_mass > 0) {
        double current_thrust = craft->thrust * craft->throttle;

        // default engine thrust direction (positive Y is "front")
        vec3 engine_thrust_direction = {0.0, 1.0, 0.0};

        // rotate thrust by spacecraft attitude
        vec3 world_thrust = quaternionRotate(craft->attitude, engine_thrust_direction);

        // apply thrust
        craft->grav_force = vec3_add(craft->grav_force, vec3_scale(world_thrust, current_thrust));
    }
}

void craft_consumeFuel(spacecraft_t* craft, const double dt) {
    if (craft->engine_on && craft->fuel_mass > 0) {
        double fuel_consumed = craft->mass_flow_rate * craft->throttle * dt;

        if (fuel_consumed > craft->fuel_mass) {
            fuel_consumed = craft->fuel_mass;
            craft->engine_on = false;
        }

        craft->fuel_mass -= fuel_consumed;
        craft->current_total_mass = craft->dry_mass + craft->fuel_mass;
    }
}

// updates the motion of the spacecraft using velocity verlet integration
void craft_updateMotion(spacecraft_t* craft, const double dt) {
    // calculate the current acceleration from the force
    craft->acc = vec3_scale(craft->grav_force, 1.0 / craft->current_total_mass);

    // update position using current velocity and acceleration
    vec3 vel_term = vec3_scale(craft->vel, dt);
    vec3 acc_term = vec3_scale(craft->acc, 0.5 * dt * dt);
    craft->pos = vec3_add(craft->pos, vec3_add(vel_term, acc_term));

    // update velocity using average of current and previous acceleration
    vec3 avg_acc = vec3_scale(vec3_add(craft->acc, craft->acc_prev), 0.5);
    craft->vel = vec3_add(craft->vel, vec3_scale(avg_acc, dt));
    craft->vel_mag = vec3_mag(craft->vel);

    // store current acceleration for next iteration
    craft->acc_prev = craft->acc;
}

// adds a spacecraft to the spacecraft array
void craft_addSpacecraft(spacecraft_properties_t* sc, const char* name,
                        vec3 pos, vec3 vel,
                        const double dry_mass, const double fuel_mass, const double thrust,
                        const double specific_impulse, const double mass_flow_rate,
                        const double attitude_angle, const double moment_of_inertia,
                        const double nozzle_gimbal_range,
                        const burn_properties_t* burns, const int num_burns) {

    // grow capacity if needed
    if (sc->count >= sc->capacity) {
        int new_capacity = sc->capacity == 0 ? 4 : sc->capacity * 2;
        spacecraft_t* temp = (spacecraft_t*)realloc(sc->spacecraft, new_capacity * sizeof(spacecraft_t));
        if (temp == NULL) {
            displayError("ERROR", "Failed to allocate memory for spacecraft");
            return;
        }
        sc->spacecraft = temp;
        sc->capacity = new_capacity;
    }

    int idx = sc->count;
    spacecraft_t* craft = &sc->spacecraft[idx];

    // allocate and copy name
    craft->name = (char*)malloc(strlen(name) + 1);
    if (craft->name == NULL) {
        displayError("ERROR", "Error: Failed to allocate memory for spacecraft name\n");
        return;
    }
#ifdef WIN32
    strcpy_s(craft->name, strlen(name) + 1, name);
#else
    strcpy(craft->name, name);
#endif

    // initialize position and velocity
    craft->pos = pos;
    craft->vel = vel;
    craft->vel_mag = vec3_mag(vel);
    craft->acc = vec3_zero();
    craft->acc_prev = vec3_zero();
    craft->grav_force = vec3_zero();

    // initialize attitude
    vec3 start_axis = {0.0, 0.0, 1.0};
    craft->attitude = quaternionFromAxisAngle(start_axis, attitude_angle);

    // initialize mass properties
    craft->dry_mass = dry_mass;
    craft->fuel_mass = fuel_mass;
    craft->current_total_mass = dry_mass + fuel_mass;

    // initialize propulsion
    craft->mass_flow_rate = mass_flow_rate;
    craft->thrust = thrust;
    craft->specific_impulse = specific_impulse;
    craft->throttle = 0.0;
    craft->engine_on = false;
    craft->nozzle_gimbal_range = nozzle_gimbal_range;
    craft->nozzle_velocity = 0.0;

    // initialize rotation/inertia
    craft->rotational_v = 0.0;
    craft->momentum = 0.0;
    craft->rotational_a = 0.0;
    craft->moment_of_inertia = moment_of_inertia;
    craft->torque = 0.0;

    // initialize SOI tracking
    craft->SOI_planet_id = 0;
    craft->SOI_planet_dist = 0;

    // initialize burn schedule
    craft->num_burns = num_burns;
    if (num_burns > 0) {
        craft->burn_properties = (burn_properties_t*)malloc(num_burns * sizeof(burn_properties_t));
        for (int i = 0; i < num_burns; i++) {
            craft->burn_properties[i] = burns[i];
        }
    } else {
        craft->burn_properties = NULL;
    }

    sc->count++;
}
