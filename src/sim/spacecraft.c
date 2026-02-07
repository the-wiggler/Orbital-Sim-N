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

            if (burn->relative_burn_target.direct) {
                // direct: attitude quaternion is pre-computed
                final_attitude = burn->burn_attitude;
            } else if (burn->relative_burn_target.absolute) {
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
        return (vec3){NAN, NAN, NAN};
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
    craft->orbital_elements.closest_r_squared = closest_r_squared;
    craft->orbital_elements.closest_planet_id = closest_planet_id;
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

vec3 craft_solveLambertProblem(const spacecraft_t* craft, const vec3 final_pos, const double time_of_flight, const body_t* central_body) {
    // position vectors relative to central body
    const vec3 r1_rel = vec3_sub(craft->pos, central_body->pos);
    const vec3 r2_rel = vec3_sub(final_pos, central_body->pos);

    const double r1_mag = vec3_mag(r1_rel);
    const double r2_mag = vec3_mag(r2_rel);

    // chord distance (distance between the two position vectors)
    const double dist = vec3_mag(vec3_sub(r2_rel, r1_rel));

    // the semi-perimeter of the triangle
    const double semi_perimeter = (r1_mag + r2_mag + dist) / 2.0;

    // transfer angle cosine
    const double cos_delta_v = vec3_dot(r1_rel, r2_rel) / (r1_mag * r2_mag);

    double a_guess = semi_perimeter / 2.0; // start at min possible value for a
    const double a_step = a_guess * 0.01; // the step by which we should increase a until we reach tof
    double tof_guess = 0.0; // calculated tof based on semi-major-axis iterations
    double prev_tof = 0.0;
    double alpha = 0.0;
    double beta = 0.0;

    // iterate through values of semi-major-axis until the calculated value of tof_guess roughly matches that of time_of_flight
    while ((tof_guess / time_of_flight) < 0.9 || (tof_guess / time_of_flight) > 1.1 ) {
        prev_tof = tof_guess;
        alpha = 2.0 * asin( sqrt( semi_perimeter / (2.0 * a_guess) ) );
        beta = 2.0 * asin( sqrt( (semi_perimeter - dist) / (2.0 * a_guess) ) );
        tof_guess = sqrt( (a_guess * a_guess * a_guess) / central_body->gravitational_parameter ) * ( alpha - sin(alpha) - (beta - sin(beta)) );
        a_guess += a_step;
        printf("tof_guess: %f | a_guess: %f\n", tof_guess, a_guess);

        // if tof is barely changing per step, it'll never converge to the target
        if (prev_tof != 0.0 && fabs(tof_guess - prev_tof) / fabs(tof_guess) < 1e-4) {
            return (vec3){INFINITY, INFINITY, INFINITY};
        }
    }

    // calculate orbital parameter
    const double sin_half_ab = sin((alpha + beta) / 2.0);
    const double orbital_parameter = ( (4.0 * a_guess) * (semi_perimeter - r2_mag) * (semi_perimeter - r1_mag) ) / (dist * dist) * (sin_half_ab * sin_half_ab);

    // calculate lagrange coefficients
    const double f_coeff = 1.0 - (r2_mag / orbital_parameter) * (1.0 - cos_delta_v);
    const double g_coeff = (r1_mag * r2_mag * sin(acos(cos_delta_v))) / sqrt(central_body->gravitational_parameter * orbital_parameter);

    // calculate v1 at departure
    const vec3 f_times_r1 = vec3_scale(r1_rel, f_coeff);
    const vec3 r2_minus_fr1 = vec3_sub(r2_rel, f_times_r1);
    const vec3 v1 = vec3_scale(r2_minus_fr1, 1.0 / g_coeff);

    // v1 is relative to central body, convert to absolute frame and compute delta-v
    const vec3 v1_absolute = vec3_add(v1, central_body->vel);
    const vec3 delta_v = vec3_sub(v1_absolute, craft->vel);

    printf("Delta V Vector: (%f, %f, %f)", delta_v.x, delta_v.y, delta_v.z); // for debug

    // returns delta v vector
    return delta_v;
}

// function that generates a burn list for the craft based on the auto target designation from the JSON file
burn_properties_t craft_createAutoTargetBurns(const sim_properties_t* sim, const int craft_id) {
    const spacecraft_t* craft = &sim->global_spacecraft.spacecraft[craft_id];
    const body_t* target_body = &sim->global_bodies.bodies[craft->auto_target_data.target_body_id];

    // TODO: grid search like method that determines optimal time of flight by solving lambert problem with like a
    // billion different random delta_t values and just picking the one with the smallest delta_v
    // one you pick a random time, then you determine where the moon will be at that time, simple 2 body problem
    // (just steal the orbital params the sim calculated for simplicity, the craft is "calculating" these)
    // compare this to previous value you calculated (probably store the delta v value from before)
    // and see if its bigger or smaller, if bigger throw it away and iterate at a different value opposite to the one
    // you just picked

    // determine the final position that the craft should be at (relative to the target body at this current point in time)
    vec3 final_pos = vec3_add((vec3){100000000, 100000000, 0}, target_body->pos);

    // estimate transfer time to set a realistic sample range
    const double r1 = vec3_mag(vec3_sub(craft->pos, target_body->pos));
    const double r2 = vec3_mag(vec3_sub(final_pos, target_body->pos));
    const double a_transfer = (r1 + r2) / 2.0;
    const double t_transfer = N_PI * sqrt((a_transfer * a_transfer * a_transfer) / target_body->gravitational_parameter);
    const double t_min = t_transfer * 0.5;
    const double t_max = t_transfer * 2.0;

    const int samples = 1000;
    double time_samples[samples];

    // populate time sample array with realistic transfer times
    const double t_step = (t_max - t_min) / (double)samples;
    time_samples[0] = t_min;
    for (int i = 1; i < samples; i++) {
        time_samples[i] = time_samples[i-1] + t_step;
    }

    // iteratively solve lambert problem with different time values to try and find the smallest delta v value
    vec3 calculated_delta_v = {INFINITY, INFINITY, INFINITY};
    for (int i = 0; i < samples; i++) {
        vec3 sample_delta_v = craft_solveLambertProblem(craft, final_pos, time_samples[i], target_body);
        if (vec3_mag_sq(sample_delta_v) < vec3_mag_sq(calculated_delta_v)) {
            calculated_delta_v = sample_delta_v;
            printf("calculated_delta_v: %f\n", vec3_mag_sq(calculated_delta_v));
        }
    }

    const double delta_v_magnitude = vec3_mag(calculated_delta_v);
    printf("Final Delta V: %f\n", delta_v_magnitude);

    // compute attitude quaternion that points the spacecraft along the delta-v direction
    // default engine thrust direction is +Y
    const vec3 default_forward = {0.0, 1.0, 0.0};
    const quaternion_t burn_attitude = quaternionFromTwoVectors(default_forward, calculated_delta_v);

    // calculate burn duration
    const double burn_duration = craft->current_total_mass * (1.0 - exp(-delta_v_magnitude / (craft->specific_impulse * STANDARD_GRAVITY))) / craft->mass_flow_rate;

    // performs the burn immediately when the function is run
    burn_properties_t burn_properties_to_achieve_target = (burn_properties_t){
        .burn_start_time = sim->window_params.sim_time,
        .burn_end_time = sim->window_params.sim_time + burn_duration,
        .throttle = 1.0,
        .burn_heading = 0.0,
        .burn_attitude = burn_attitude,
        .burn_target_id = craft->auto_target_data.target_body_id,
        .auto_burn = true,
        .auto_burn_final_pos = final_pos,
        .relative_burn_target = { .direct = true }
    };

    return burn_properties_to_achieve_target;
}

// adds a spacecraft to the spacecraft array
// NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
void craft_addSpacecraft(spacecraft_properties_t* global_spacecraft, const char* name, const vec3 pos, const vec3 vel, const double dry_mass, const double fuel_mass, const double thrust, const double specific_impulse, const double mass_flow_rate, const double attitude_angle, const double moment_of_inertia, const double nozzle_gimbal_range, const double nozzle_velocity, const burn_properties_t* burns, const int num_burns) {

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
    craft->nozzle_velocity = nozzle_velocity;

    // initialize rotation/inertia
    craft->rotational_v = 0.0;
    craft->momentum = 0.0;
    craft->rotational_a = 0.0;
    craft->moment_of_inertia = moment_of_inertia;
    craft->torque = 0.0;

    // initialize SOI tracking
    craft->orbital_elements.SOI_planet_id = 0;
    craft->orbital_elements.closest_r_squared = INFINITY;
    craft->orbital_elements.closest_planet_id = 0;

    craft->orbital_elements.apoapsis = 0.0;
    craft->orbital_elements.periapsis = 0.0;
    craft->orbital_elements.semi_major_axis = 0.0;
    craft->orbital_elements.eccentricity = 0.0;

    // initialize burn schedule
    if (num_burns > MAX_BURNS_PER_SPACECRAFT) {
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
