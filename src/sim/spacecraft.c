#include "spacecraft.h"
#include "../globals.h"
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "../types.h"

#include "../math/matrix.h"

void displayError(const char* title, const char* message);

// check and activate burns
// NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
void craft_checkBurnSchedule(spacecraft_t* craft, const body_properties_t* global_bodies, const double sim_time) {
    bool burn_active = false;
    for (int j = 0; j < craft->num_burns; j++) {
        const burn_properties_t* burn = &craft->burn_properties[j];
        // check if within the burn window
        if (sim_time >= burn->burn_start_time && sim_time < burn->burn_end_time && craft->fuel_mass > 0) {
            craft->engine_on = true;
            craft->throttle = burn->throttle;

            // calculate attitude based on burn type
            quaternion_t final_attitude = {0};
            const int target_id = burn->burn_target_id;
            const body_t* target = &global_bodies->bodies[target_id];

            if (burn->relative_burn_target.absolute) {
                // absolute: heading is in absolute space coordinates
                const vec3 absolute_axis = {0.0, 0.0, 1.0};
                final_attitude = quaternionFromAxisAngle(absolute_axis, burn->burn_heading);
            } else if (burn->relative_burn_target.tangent) {
                // tangent: heading is relative to the velocity vector
                const vec3 rel_vel = vec3_sub(craft->vel, target->vel);
                const vec3 default_forward = {0.0, 1.0, 0.0};
                const quaternion_t base_rotation = quaternionFromTwoVectors(default_forward, rel_vel);

                if (burn->burn_heading != 0.0) {
                    const vec3 rotation_axis = vec3_normalize(rel_vel);
                    const quaternion_t offset_rotation = quaternionFromAxisAngle(rotation_axis, burn->burn_heading);
                    final_attitude = quaternionMul(base_rotation, offset_rotation);
                } else {
                    final_attitude = base_rotation;
                }
            } else if (burn->relative_burn_target.normal) {
                // normal: heading is perpendicular to the orbital plane
                const vec3 rel_pos = vec3_sub(craft->pos, target->pos);
                const vec3 rel_vel = vec3_sub(craft->vel, target->vel);
                const vec3 normal_direction = cross_product_vec3(rel_pos, rel_vel);
                const vec3 default_forward = {0.0, 1.0, 0.0};
                const quaternion_t base_rotation = quaternionFromTwoVectors(default_forward, normal_direction);

                if (burn->burn_heading != 0.0) {
                    const vec3 rotation_axis = vec3_normalize(normal_direction);
                    const quaternion_t offset_rotation = quaternionFromAxisAngle(rotation_axis, burn->burn_heading);
                    final_attitude = quaternionMul(base_rotation, offset_rotation);
                } else {
                    final_attitude = base_rotation;
                }
            } else {
                displayError("ERROR", "Failed at determining burn type. If you see this the dev sucks at coding lol");
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
// NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
vec3 craft_calculateGravForce(const sim_properties_t* sim, const int craft_idx, const int body_idx) {
    const spacecraft_t* craft = &sim->global_spacecraft.spacecraft[craft_idx];
    const body_t* body = &sim->global_bodies.bodies[body_idx];

    // calculate the distance between the spacecraft and the body
    const vec3 delta_pos = vec3_sub(body->pos, craft->pos);
    const double r_squared = vec3_mag_sq(delta_pos);
    const double radius = sqrt(r_squared);

    // planet collision logic
    if (radius < body->radius) {
        char err_txt[MAX_ERR_SIZE];
        snprintf(err_txt, sizeof(err_txt), "Warning: %s has collided with %s\n\nResetting Simulation...", craft->name, body->name);
        displayError("PLANET COLLISION", err_txt);
        return vec3_zero();
    }

    // force = (G * m1 * m2) * delta / r^3
    const double r_cubed = r_squared * radius;
    const double force_factor = (G_C * craft->current_total_mass * body->mass) / r_cubed;

    // apply the force to the craft
    const vec3 force = vec3_scale(delta_pos, force_factor);

    return force;
}

// calculates the J2 perturbation force applied to a spacecraft
// accounts for the body's axial tilt by transforming to body-fixed frame
// derived from this equation: https://www.vcalc.com/equation/?uuid=1e5aa6ea-95a3-11e7-9770-bc764e2038f2 -- sorry this isn't a reliable source, but it's the best representation of the equation that I could easily find
// NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
vec3 craft_calculateJ2Force(const sim_properties_t* sim, const int craft_idx, const int body_idx) {
    const spacecraft_t* craft = &sim->global_spacecraft.spacecraft[craft_idx];
    const body_t* body = &sim->global_bodies.bodies[body_idx];

    // position relative to body in world space
    const vec3 rel_pos_world = vec3_sub(craft->pos, body->pos);

    // transform to body-fixed frame where Z is aligned with spin axis
    // use conjugate of attitude to go from world -> body frame
    const quaternion_t q_inv = quaternionConjugate(body->attitude);
    const vec3 rel_pos_body = quaternionRotate(q_inv, rel_pos_world);

    // apply J2 formula in body frame of reference
    const double radius = vec3_mag(rel_pos_body);
    const double r_squared = radius * radius;
    const double r_5 = r_squared * r_squared * radius;

    const double term1 = craft->current_total_mass * (body->J2_per_coeff_numerator / r_5);
    const double z_r_ratio_sq = (rel_pos_body.z * rel_pos_body.z) / r_squared;
    const double term2 = 5.0 * z_r_ratio_sq;

    // J2 force in body frame
    vec3 force_body;
    force_body.x = term1 * (term2 - 1.0) * rel_pos_body.x;
    force_body.y = term1 * (term2 - 1.0) * rel_pos_body.y;
    force_body.z = term1 * (term2 - 3.0) * rel_pos_body.z;

    // transform force back to world frame
    return quaternionRotate(body->attitude, force_body);
}

// updates the ID and distance of the closest planet
// (this should be used when initially spawning a craft because the grav calculations do this exact calculation by default)
// this is probably executed when the JSON is loaded
void craft_findClosestPlanet(spacecraft_t* craft, const body_properties_t* global_bodies) {
    double closest_r_squared = INFINITY;
    int closest_planet_id = 0;
    for (int i = 0; i < global_bodies->count; i++) {
        // calculate the distance between the spacecraft and the body
        const vec3 delta_pos = vec3_sub(global_bodies->bodies[i].pos, craft->pos);
        const double r_squared = vec3_mag_sq(delta_pos);
        if (r_squared < closest_r_squared) {
            closest_r_squared = r_squared;
            closest_planet_id = i;
        }
    }
    craft->oe.closest_r_squared = closest_r_squared;
    craft->oe.closest_planet_id = closest_planet_id;
}

// applies thrust force based on current attitude
void craft_applyThrust(spacecraft_t* craft) {
    if (craft->engine_on && craft->fuel_mass > 0) {
        const double current_thrust = craft->thrust * craft->throttle;

        // default engine thrust direction (positive Y is "front")
        const vec3 engine_thrust_direction = {0.0, 1.0, 0.0};

        // rotate thrust by spacecraft attitude
        const vec3 world_thrust = quaternionRotate(craft->attitude, engine_thrust_direction);

        // apply thrust
        craft->grav_force = vec3_add(craft->grav_force, vec3_scale(world_thrust, current_thrust));
    }
}

void craft_consumeFuel(spacecraft_t* craft, const double delta_t) {
    if (craft->engine_on && craft->fuel_mass > 0) {
        double fuel_consumed = craft->mass_flow_rate * craft->throttle * delta_t;

        if (fuel_consumed > craft->fuel_mass) {
            fuel_consumed = craft->fuel_mass;
            craft->engine_on = false;//
        }

        craft->fuel_mass -= fuel_consumed;
        craft->current_total_mass = craft->dry_mass + craft->fuel_mass;
    }
}

// adds a spacecraft to the spacecraft array
// NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
void craft_addSpacecraft(spacecraft_properties_t* global_spacecraft, const char* name, const vec3 pos, const vec3 vel, const double dry_mass, const double fuel_mass, const double thrust, const double specific_impulse, const double mass_flow_rate, const double attitude_angle, const double moment_of_inertia, const double nozzle_gimbal_range, const burn_properties_t* burns, const int num_burns) {

    // check bounds
    if (global_spacecraft->count >= MAX_SPACECRAFT) {
        char err_txt[MAX_ERR_SIZE];
        snprintf(err_txt, sizeof(err_txt), "Cannot add spacecraft '%s': Maximum of %d spacecraft reached", name, MAX_SPACECRAFT);
        displayError("ERROR", err_txt);
        return;
    }

    const int idx = global_spacecraft->count;
    spacecraft_t* craft = &global_spacecraft->spacecraft[idx];

    // validate and copy name
    const size_t name_len = strlen(name);
    if (name_len >= MAX_NAME_LENGTH) {
        char err_txt[MAX_ERR_SIZE];
        snprintf(err_txt, sizeof(err_txt), "Warning: Spacecraft name '%s' truncated to %d characters", name, MAX_NAME_LENGTH - 1);
        displayError("WARNING", err_txt);
    }
#ifdef _WIN32
    strncpy_s(craft->name, MAX_NAME_LENGTH, name, MAX_NAME_LENGTH - 1);
#else
    strncpy(craft->name, name, MAX_NAME_LENGTH - 1);
#endif
    craft->name[MAX_NAME_LENGTH - 1] = '\0';  // ensure null termination

    // initialize position and velocity
    craft->pos = pos;
    craft->vel = vel;
    craft->vel_mag = vec3_mag(vel);
    craft->acc = vec3_zero();
    craft->acc_prev = vec3_zero();
    craft->grav_force = vec3_zero();

    // initialize attitude
    const vec3 start_axis = {0.0, 0.0, 1.0};
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
    craft->oe.SOI_planet_id = 0;
    craft->oe.closest_r_squared = INFINITY;
    craft->oe.closest_planet_id = 0;

    craft->oe.apoapsis = 0.0;
    craft->oe.periapsis = 0.0;
    craft->oe.semi_major_axis = 0.0;
    craft->oe.eccentricity = 0.0;

    // initialize burn schedule
    if (num_burns > MAX_BURNS_PER_SPACECRAFT) {
        char err_txt[MAX_ERR_SIZE];
        snprintf(err_txt, sizeof(err_txt), "Warning: Spacecraft '%s' has %d burns, exceeding maximum of %d. Truncating.",
                 name, num_burns, MAX_BURNS_PER_SPACECRAFT);
        displayError("WARNING", err_txt);
        craft->num_burns = MAX_BURNS_PER_SPACECRAFT;
    } else {
        craft->num_burns = num_burns;
    }

    // copy burns into fixed array
    for (int i = 0; i < craft->num_burns; i++) {
        craft->burn_properties[i] = burns[i];
    }

    global_spacecraft->count++;
}
