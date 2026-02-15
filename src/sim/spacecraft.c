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

// solves the lambert problem with position inputs and delta v as an output to get to the desired position
// NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
vec3 craft_solveLambertProblem(const vec3 craft_pos, const vec3 craft_vel, const vec3 final_pos, const double time_of_flight, const body_t* central_body) {
    // position vectors relative to central body
    const vec3 r1_rel = vec3_sub(craft_pos, central_body->pos);
    const vec3 r2_rel = vec3_sub(final_pos, central_body->pos);

    const double r1_mag = vec3_mag(r1_rel);
    const double r2_mag = vec3_mag(r2_rel);

    // chord distance (distance between the two position vectors)
    const double dist = vec3_mag(vec3_sub(r2_rel, r1_rel));

    // the semi-perimeter of the triangle
    const double semi_perimeter = (r1_mag + r2_mag + dist) / 2.0;

    // transfer angle cosine
    const double cos_delta_v = vec3_dot(r1_rel, r2_rel) / (r1_mag * r2_mag);

    // determine transfer direction (short-way vs long-way) for prograde transfers
    const vec3 transfer_cross = cross_product_vec3(r1_rel, r2_rel);
    const bool long_way = transfer_cross.z < 0.0;

    double a_guess = semi_perimeter / 2.0; // start at min possible value for a
    const double a_step = a_guess * 0.01; // the step by which we should increase a until we reach tof
    double tof_guess = 0.0; // calculated tof based on semi-major-axis iterations
    double prev_tof = 0.0;
    double alpha = 0.0;
    double beta = 0.0;

    // iterate through values of semi-major-axis until the calculated value of tof_guess closely matches that of time_of_flight
    while ((tof_guess / time_of_flight) < 0.99999 || (tof_guess / time_of_flight) > 1.00001 ) {
        prev_tof = tof_guess;
        alpha = 2.0 * asin( sqrt( semi_perimeter / (2.0 * a_guess) ) );
        beta = 2.0 * asin( sqrt( (semi_perimeter - dist) / (2.0 * a_guess) ) );
        if (long_way) { beta = -beta; }
        tof_guess = sqrt( (a_guess * a_guess * a_guess) / central_body->gravitational_parameter ) * ( alpha - sin(alpha) - (beta - sin(beta)) );
        //printf("tof_guess: %f | a_guess: %f\n", tof_guess, a_guess); // debug info -- only turn this on if necessary it prints A LOT of lines

        // if tof is barely changing per step, it'll never converge to the target
        if (prev_tof != 0.0 && fabs(tof_guess - prev_tof) / fabs(tof_guess) < 1e-4) {
            return (vec3){INFINITY, INFINITY, INFINITY};
        }
        a_guess += a_step;
    }
    a_guess -= a_step; // correct for the extra increment after convergence

    // calculate orbital parameter
    const double sin_half_ab = sin((alpha + beta) / 2.0);
    const double orbital_parameter = ( (4.0 * a_guess) * (semi_perimeter - r2_mag) * (semi_perimeter - r1_mag) ) / (dist * dist) * (sin_half_ab * sin_half_ab);

    // calculate lagrange coefficients
    const double f_coeff = 1.0 - ((r2_mag / orbital_parameter) * (1.0 - cos_delta_v));
    const double sin_delta_nu = long_way ? -sqrt(1.0 - (cos_delta_v * cos_delta_v)) : sqrt(1.0 - (cos_delta_v * cos_delta_v));
    const double g_coeff = (r1_mag * r2_mag * sin_delta_nu) / sqrt(central_body->gravitational_parameter * orbital_parameter);

    // calculate vel1 at departure
    const vec3 f_times_r1 = vec3_scale(r1_rel, f_coeff);
    const vec3 r2_minus_fr1 = vec3_sub(r2_rel, f_times_r1);
    const vec3 vel1 = vec3_scale(r2_minus_fr1, 1.0 / g_coeff);

    // vel1 is relative to central body, convert to absolute frame and compute delta-v
    const vec3 v1_absolute = vec3_add(vel1, central_body->vel);
    const vec3 delta_v = vec3_sub(v1_absolute, craft_vel);

    // returns delta v vector
    return delta_v;
}

// determines the time for a hohmann transfer
// NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
double craft_calcualteHohmannTransTime(const double craft_orbit_radius, const double target_orbit_radius, const double central_grav_param) {
    const double a_transfer = (craft_orbit_radius + target_orbit_radius) / 2.0;
    const double transfer_time = N_PI * sqrt((a_transfer * a_transfer * a_transfer) / central_grav_param);
    return transfer_time;
}

// solves kepler's equation M = E - e*sin(E) for eccentric anomaly
static double solveKeplerEquation(double mean_anomaly, const double eccentricity) {
    mean_anomaly = fmod(mean_anomaly, TWO_PI);
    if (mean_anomaly < 0.0) { mean_anomaly += TWO_PI; }

    double eccentric_anomaly = mean_anomaly;
    for (int iter = 0; iter < 50; iter++) {
        const double d_eccentric_anomaly = (eccentric_anomaly - eccentricity * sin(eccentric_anomaly) - mean_anomaly) / (1.0 - eccentricity * cos(eccentric_anomaly));
        eccentric_anomaly -= d_eccentric_anomaly;
        if (fabs(d_eccentric_anomaly) < 1e-12) { break; }
    }
    return eccentric_anomaly;
}

// propagates an orbit forward by dt seconds using kepler's equation and returns position and velocity
// in the absolute ref frame. central_pos/vel are the central body's current position/velocity
// thank you clanker 😍
// NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
static void propagateOrbitState(const orbital_elements_t* orbital_elements, const double grav_param, const vec3 central_pos, const vec3 central_vel, const double delta_t, vec3* out_pos, vec3* out_vel) {
    const double sma = orbital_elements->semi_major_axis;
    const double eccentricity = orbital_elements->eccentricity;
    const double inc = orbital_elements->inclination;
    const double Omega = orbital_elements->ascending_node;
    const double omega = orbital_elements->arg_periapsis;
    const double nu0 = orbital_elements->true_anomaly;
    const double mean_motion = sqrt(grav_param / (sma * sma * sma));

    const double initial_E = atan2(sqrt(1.0 - (eccentricity * eccentricity)) * sin(nu0), eccentricity + cos(nu0));
    const double initial_M = initial_E - (eccentricity * sin(initial_E));

    const double eccentric_anomaly = solveKeplerEquation(initial_M + (mean_motion * delta_t), eccentricity);

    const double true_anomaly = atan2(sqrt(1.0 - (eccentricity * eccentricity)) * sin(eccentric_anomaly), cos(eccentric_anomaly) - eccentricity);
    const double radius = sma * (1.0 - (eccentricity * cos(eccentric_anomaly)));

    const vec3 z_axis = {0.0, 0.0, 1.0};
    const vec3 x_axis = {1.0, 0.0, 0.0};
    const quaternion_t q_frame = quaternionMul(quaternionMul(
        quaternionFromAxisAngle(z_axis, Omega),
        quaternionFromAxisAngle(x_axis, inc)),
        quaternionFromAxisAngle(z_axis, omega));

    const vec3 pos_orb = {radius * cos(true_anomaly), radius * sin(true_anomaly), 0.0};
    const vec3 pos_rel = quaternionRotate(q_frame, pos_orb);
    *out_pos = vec3_add(pos_rel, central_pos);

    const double semi_latus = sma * (1.0 - eccentricity * eccentricity);
    const double sqrt_mu_over_p = sqrt(grav_param / semi_latus);
    const vec3 vel_orb = {-sqrt_mu_over_p * sin(true_anomaly), sqrt_mu_over_p * (eccentricity + cos(true_anomaly)), 0.0};
    const vec3 vel_rel = quaternionRotate(q_frame, vel_orb);
    *out_vel = vec3_add(vel_rel, central_vel);
}

// function that generates a burn list for the craft based on the optimal delta v for the final desired position.
// IMPORTANT: THIS FUNCTION SHOULD ONLY BE RUN ON THE MAIN THREAD, NOT THE SIM THREAD. THE RESULT IS UNKNOWN AND POTENTIALLY DANGEROUS IF RUN ON THE SIM THREAD
burn_properties_t craft_autoDeltaVOptimization(const sim_properties_t* sim, const int craft_id) {
    const spacecraft_t* craft = &sim->global_spacecraft.spacecraft[craft_id];
    const body_t* central_body = &sim->global_bodies.bodies[craft->orbital_elements.SOI_planet_id];
    const body_t* target_body = &sim->global_bodies.bodies[craft->target_body_id];

    printf("Started Burn Optimization...\n");

    // grid search like method that determines optimal delta v by solving lambert problem with like a
    // billion different random delta_t values and just picking the one with the smallest delta_v
    const int ORBIT_INCREMENT_RES = 500;
    const double grav_param = central_body->gravitational_parameter;

    // determine the orbital period for the craft
    const double craft_orbital_period = 2.0 * N_PI * sqrt((craft->orbital_elements.semi_major_axis *
        craft->orbital_elements.semi_major_axis * craft->orbital_elements.semi_major_axis) / grav_param);
    const double craft_orbit_time_step = craft_orbital_period / ORBIT_INCREMENT_RES;

    // propagate craft orbit using kepler's equation to get position AND velocity at each time step
    vec3 craft_locations[ORBIT_INCREMENT_RES];
    vec3 craft_velocities[ORBIT_INCREMENT_RES];

    for (int i = 0; i < ORBIT_INCREMENT_RES; i++) {
        const double delta_t = i * craft_orbit_time_step;
        propagateOrbitState(&craft->orbital_elements, grav_param, central_body->pos, central_body->vel,
                           delta_t, &craft_locations[i], &craft_velocities[i]);
    }

    // target body orbital radius (for hohmann estimate)
    const double target_orbit_radius = vec3_mag(vec3_sub(target_body->pos, central_body->pos));

    // increment over different times in the craft's orbit to determine where in the orbit the most optimal delta v exists
    vec3 best_delta_v = {INFINITY, INFINITY, INFINITY};
    double expected_time_of_contact = 0;
    int best_i = 0;
    for (int i = 0; i < ORBIT_INCREMENT_RES; i++) {
        const double dist_from_central_body = vec3_mag(vec3_sub(craft_locations[i], central_body->pos));

        // 1- determine ToF bounds of grid search using hohmann transfer
        const double hohmann_time = craft_calcualteHohmannTransTime(dist_from_central_body, target_orbit_radius, grav_param);
        const double tof_lower = hohmann_time * 0.8;
        const double tof_upper = hohmann_time * 1.2;

        // 2- do a rough grid search at current point in the orbit to find where the planet is and lowest delta v
        const int SAMPLES = 500;
        const double t_step = (tof_upper - tof_lower) / (double)SAMPLES;

        // populate time sample array with potential realistic transfer times
        // and determine where the planet's position is at every potential transfer time
        double time_samples[SAMPLES];
        time_samples[0] = tof_lower;
        vec3 target_body_location[SAMPLES];

        for (int j = 0; j < SAMPLES; j++) {
            if (j > 0) {time_samples[j] = time_samples[j - 1] + t_step;}

            // propagate target body forward by time for craft to reach orbit point i + transfer time
            const double total_dt = (i * craft_orbit_time_step) + time_samples[j];
            vec3 tgt_vel_unused;
            propagateOrbitState(&target_body->oe, grav_param, central_body->pos, central_body->vel,
                               total_dt, &target_body_location[j], &tgt_vel_unused);
        }

        // determine how much delta v is needed for this circumstance//
        vec3 calculated_delta_v = {INFINITY, INFINITY, INFINITY};
        double best_tof_at_i = 0.0;
        for (int j = 0; j < SAMPLES; j++) {
            // delta v to get to that planets location at given time by solving lambert problem
            vec3 sample_delta_v = craft_solveLambertProblem(craft_locations[i], craft_velocities[i], target_body_location[j], time_samples[j], central_body);

            // if this sample delta v is smaller than our smallest yet, we make this one the new smallest ever (you are the smallest delta_v ever!)
            if (vec3_mag_sq(sample_delta_v) < vec3_mag_sq(calculated_delta_v)) {
                calculated_delta_v = sample_delta_v;
                best_tof_at_i = time_samples[j];
            }
        }

        // 3- check if this calculated optimal delta v is smaller than the one previously calculated at a different point in the orbit
        if (vec3_mag_sq(calculated_delta_v) < vec3_mag_sq(best_delta_v)) {
            best_delta_v = calculated_delta_v;
            best_i = i;
            expected_time_of_contact = best_tof_at_i;
        }
    
        // advance to a farther time in the orbit to see if any of the possible delta v values are better there, because the ideal time
        // to burn COULD be in the future and not immediately
        printf("\rProgress: %.1f%%", (float)i / (float)ORBIT_INCREMENT_RES * 100.0F);
        fflush(stdout);
    }
    printf("\n");

    const double best_delta_v_magnitude = vec3_mag(best_delta_v);

    // compute attitude quaternion that points the spacecraft along the delta-v direction
    // default engine thrust direction is +Y
    const vec3 default_forward = {0.0, 1.0, 0.0};
    const quaternion_t burn_attitude = quaternionFromTwoVectors(default_forward, best_delta_v);

    // calculate burn duration
    const double burn_duration = craft->current_total_mass * (1.0 - exp(-best_delta_v_magnitude / (craft->specific_impulse * STANDARD_GRAVITY))) / craft->mass_flow_rate;

    // schedule the burn centered on the optimal point in the orbit
    const double burn_start_time = sim->window_params.sim_time + (best_i * craft_orbit_time_step) - (burn_duration / 2.0);

    printf("Final Delta V: %f m/s | Burn occurs in %e s | Transfer time: %e s | Burn time: %f\n",
       best_delta_v_magnitude, best_i * craft_orbit_time_step, expected_time_of_contact, burn_duration);

    burn_properties_t burn_properties_to_achieve_target = (burn_properties_t){
        .burn_start_time = burn_start_time,
        .burn_end_time = burn_start_time + burn_duration,
        .throttle = 1.0,
        .burn_heading = 0.0,
        .burn_attitude = burn_attitude,
        .burn_target_id = craft->target_body_id,
        .auto_burn = true,
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
