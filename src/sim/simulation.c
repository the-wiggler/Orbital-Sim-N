#include "simulation.h"
#include "../globals.h"
#include "../sim/bodies.h"
#include "../sim/spacecraft.h"
#include "../math/matrix.h"
#include <math.h>
#include "../types.h"


// calculate total system energy for all bodies
double calculateTotalSystemEnergy(const sim_properties_t* sim) {
    const body_properties_t* global_bodies = &sim->global_bodies;
    const spacecraft_properties_t* global_spacecraft = &sim->global_spacecraft;

    double total_kinetic = 0.0;
    double total_potential = 0.0;

    // calculate kinetic energy for all bodies
    for (int i = 0; i < global_bodies->count; i++) {
        const body_t* body = &global_bodies->bodies[i];

        total_kinetic += 0.5 * body->mass * body->vel_mag * body->vel_mag;
    }

    // calculate kinetic energy for all spacecraft
    for (int i = 0; i < global_spacecraft->count; i++) {
        const spacecraft_t* craft = &global_spacecraft->spacecraft[i];

        total_kinetic += 0.5 * craft->current_total_mass * craft->vel_mag * craft->vel_mag;
    }

    // calculate potential energy between all body pairs
    for (int i = 0; i < global_bodies->count; i++) {
        for (int j = i + 1; j < global_bodies->count; j++) {
            const body_t* body_i = &global_bodies->bodies[i];
            const body_t* body_j = &global_bodies->bodies[j];
            const vec3 delta = vec3_sub(body_j->pos, body_i->pos);
            const double radius = vec3_mag(delta);
            if (radius > 0) {
                total_potential += -(G_C * body_i->mass * body_j->mass) / radius;
            }
        }
    }

    // calculate potential energy between spacecraft and bodies
    for (int i = 0; i < global_spacecraft->count; i++) {
        for (int j = 0; j < global_bodies->count; j++) {
            const spacecraft_t* craft = &global_spacecraft->spacecraft[i];
            const body_t* body = &global_bodies->bodies[j];
            const vec3 delta = vec3_sub(body->pos, craft->pos);
            const double radius = vec3_mag(delta);
            if (radius > 0) {
                total_potential += -(G_C * craft->current_total_mass * body->mass) / radius;
            }
        }
    }

    return total_kinetic + total_potential;
}

// calculates orbital elements (this probably needs to be optimized somehow at some point because this seems very resource heavy)
// NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
void calculateOrbitalElements(orbital_elements_t* target_orbital_elements, const vec3* target_pos, const vec3* target_vel, const body_t* orbiting_body) {
    // first, the initial properties of the craft relative to the target planet should be calculated
    const vec3 c_pos     = vec3_sub(*target_pos, orbiting_body->pos); // position vector
    const vec3 c_vel     = vec3_sub(*target_vel, orbiting_body->vel); // velocity vector
    const double c_r     = vec3_mag(c_pos); // distance
    const double c_speed = vec3_mag(c_vel);
    const double grav_p  = G_C * orbiting_body->mass; // gravitational parameter
    const vec3 c_h       = vec3_cross(c_pos, c_vel); // specific angular momentum
    const vec3 z_vec     = { 0, 0, 1 };
    const vec3 c_n       = vec3_cross(z_vec, c_h); // ascending node vector

    const vec3 term1     = vec3_scalar_div(vec3_cross(c_vel, c_h), grav_p); // first term of e_vec
    const vec3 term2     = vec3_scalar_div(c_pos, c_r); // second term of e_vec
    const vec3 e_vec     = vec3_sub(term1, term2); // eccentricity vector
    target_orbital_elements->specific_E    = ((c_speed * c_speed) / 2) - (grav_p / c_r); // specific orbital energy

    // orbital elements
    target_orbital_elements->semi_major_axis = -1.0 * (grav_p / (2 * target_orbital_elements->specific_E));
    target_orbital_elements->eccentricity = vec3_mag(e_vec);

    const double h_mag = vec3_mag(c_h);
    target_orbital_elements->inclination = acos(c_h.z / h_mag); // the angle between the orbital and equatorial planes

    // longitude of ascending node -- the angle from the vernal equinox vector to the ascending node on the equatorial plane
    const double n_mag = vec3_mag(c_n);
    if (n_mag > NEAR_ZERO) {
        target_orbital_elements->ascending_node = atan2(c_n.y, c_n.x);
        if (target_orbital_elements->ascending_node < 0) {
            target_orbital_elements->ascending_node += (2 * N_PI);
        }
    } else {
        target_orbital_elements->ascending_node = 0.0; // undefined for equatorial orbits (probably unlikely to happen perfectly)
    }

    // argument of periapsis -- the angle measured between the ascending node and the perigee
    if (target_orbital_elements->eccentricity > NEAR_ZERO && n_mag > NEAR_ZERO) {
        const double cos_omega = vec3_dot(c_n, e_vec) / (n_mag * target_orbital_elements->eccentricity);
        target_orbital_elements->arg_periapsis = acos(fmax(-1.0, fmin(1.0, cos_omega)));
        if (e_vec.z < 0) {
            target_orbital_elements->arg_periapsis = (2 * N_PI) - target_orbital_elements->arg_periapsis;
        }
    } else if (target_orbital_elements->eccentricity > NEAR_ZERO) {
        // equatorial orbit, use longitude of periapsis
        target_orbital_elements->arg_periapsis = atan2(e_vec.y, e_vec.x);
        if (target_orbital_elements->arg_periapsis < 0) {
            target_orbital_elements->arg_periapsis += (2 * N_PI);
        }
    } else {
        target_orbital_elements->arg_periapsis = 0.0; // undefined for circular orbits
    }

    // true anomaly -- the angle between perigee and satellite in the orbital plane at a specific time
    if (target_orbital_elements->eccentricity > NEAR_ZERO) {
        const double cos_nu = vec3_dot(e_vec, c_pos) / (target_orbital_elements->eccentricity * c_r);
        target_orbital_elements->true_anomaly = acos(fmax(-1.0, fmin(1.0, cos_nu)));
        if (vec3_dot(c_pos, c_vel) < 0) {
            target_orbital_elements->true_anomaly = (2 * N_PI) - target_orbital_elements->true_anomaly;
        }
    } else {
        // circular orbit, use argument of latitude
        if (n_mag > NEAR_ZERO) {
            const double cos_u = vec3_dot(c_n, c_pos) / (n_mag * c_r);
            target_orbital_elements->true_anomaly = acos(fmax(-1.0, fmin(1.0, cos_u)));
            if (c_pos.z < 0) {
                target_orbital_elements->true_anomaly = (2 * N_PI) - target_orbital_elements->true_anomaly;
            }
        } else {
            // equatorial and circular, use true longitude
            target_orbital_elements->true_anomaly = atan2(c_pos.y, c_pos.x);
            if (target_orbital_elements->true_anomaly < 0) {
                target_orbital_elements->true_anomaly += (2 * N_PI);
            }
        }
    }
}

// update orbital elements and SOI tracking for bodies and spacecraft in the system
// this is separated from physics calculations and can be called less frequently for performance
// NOLINTNEXTLINE(readability-function-cognitive-complexity)
void updateSystemOrbitalElements(sim_properties_t* sim) {
    body_properties_t* global_bodies = &sim->global_bodies;
    spacecraft_properties_t* global_spacecraft = &sim->global_spacecraft;

    // find the closest planet and update SOI for all bodies
    for (int i = 0; i < global_bodies->count; i++) {
        body_t* body = &global_bodies->bodies[i];
        body->oe.closest_r_squared = INFINITY;
        body->oe.closest_planet_id = 0;

        // find the closest body to this one
        for (int j = 0; j < global_bodies->count; j++) {
            if (i == j) { continue; }

            const body_t* other_body = &global_bodies->bodies[j];
            const vec3 delta = vec3_sub(other_body->pos, body->pos);
            const double r_squared = vec3_mag_sq(delta);

            if (r_squared < body->oe.closest_r_squared) {
                body->oe.closest_r_squared = r_squared;
                body->oe.closest_planet_id = j;

                // check if within SOI
                const double radius = sqrt(r_squared);
                if (radius <= other_body->SOI_radius) {
                    body->oe.SOI_planet_id = j;
                }
            }
        }
    }

    // find the closest planet and update SOI for all spacecraft
    for (int i = 0; i < global_spacecraft->count; i++) {
        spacecraft_t* craft = &global_spacecraft->spacecraft[i];
        craft->orbital_elements.closest_r_squared = INFINITY;
        craft->orbital_elements.closest_planet_id = 0;

        // find the closest body to this spacecraft
        for (int j = 0; j < global_bodies->count; j++) {
            const body_t* body = &global_bodies->bodies[j];
            const vec3 delta = vec3_sub(body->pos, craft->pos);
            const double r_squared = vec3_mag_sq(delta);

            if (r_squared < craft->orbital_elements.closest_r_squared) {
                craft->orbital_elements.closest_r_squared = r_squared;
                craft->orbital_elements.closest_planet_id = j;

                // check if within SOI
                const double radius = sqrt(r_squared);
                if (radius <= body->SOI_radius) {
                    craft->orbital_elements.SOI_planet_id = j;
                }
            }
        }
    }

    // check if bodies have exited their current SOI
    for (int i = 0; i < global_bodies->count; i++) {
        body_t* body = &global_bodies->bodies[i];
        if (body->oe.SOI_planet_id > 0 && body->oe.SOI_planet_id < global_bodies->count) {
            const body_t* soi_body = &global_bodies->bodies[body->oe.SOI_planet_id];
            const vec3 delta = vec3_sub(soi_body->pos, body->pos);
            const double dist = vec3_mag(delta);
            if (dist > soi_body->SOI_radius) {
                // exited SOI, fall back to the closest body
                body->oe.SOI_planet_id = body->oe.closest_planet_id;
            }
        }
    }

    // calculate orbital elements for bodies (skip central body at index 0)
    for (int i = 1; i < global_bodies->count; i++) {
        body_t* body = &global_bodies->bodies[i];
        if (body->oe.SOI_planet_id != i && body->oe.SOI_planet_id >= 0 && body->oe.SOI_planet_id < global_bodies->count) {
            calculateOrbitalElements(&body->oe, &body->pos, &body->vel, &global_bodies->bodies[body->oe.SOI_planet_id]);
        }
    }

    // check if spacecraft have exited their current SOI and calculate orbital elements
    for (int i = 0; i < global_spacecraft->count; i++) {
        spacecraft_t* craft = &global_spacecraft->spacecraft[i];

        // check if spacecraft has exited its current SOI
        if (craft->orbital_elements.SOI_planet_id > 0 && craft->orbital_elements.SOI_planet_id < global_bodies->count) {
            const body_t* soi_body = &global_bodies->bodies[craft->orbital_elements.SOI_planet_id];
            const vec3 delta = vec3_sub(soi_body->pos, craft->pos);
            const double dist = vec3_mag(delta);
            if (dist > soi_body->SOI_radius) {
                // exited SOI, fall back to the closest body
                craft->orbital_elements.SOI_planet_id = craft->orbital_elements.closest_planet_id;
            }
        }
    }
}

// reset the simulation by removing all bodies from the system
void resetSim(sim_properties_t* sim) {
    body_properties_t* global_bodies = &sim->global_bodies;
    spacecraft_properties_t* global_spacecraft = &sim->global_spacecraft;
    window_params_t* window_params = &sim->window_params;

    // reset simulation time
    window_params->sim_time = 0;
    window_params->reset_sim = false;

    // reset all bodies
    global_bodies->count = 0;

    // reset all spacecraft
    global_spacecraft->count = 0;
}




// probably the most important function in the entire program. This is responsible for updating the position of literally
// every single planet each time step. It uses a method called Velocity Verlet Integration. I apologize for the function
// being so complex; it made sense to me for it all to be in one function for ease of debugging/understanding the math
// https://en.wikipedia.org/wiki/Verlet_integration#Velocity_Verlet
// NOLINTNEXTLINE(readability-function-cognitive-complexity)
void runCalculations(sim_properties_t* sim) {
    body_properties_t* global_bodies = &sim->global_bodies;
    spacecraft_properties_t* global_spacecraft = &sim->global_spacecraft;
    window_params_t* window_params = &sim->window_params;

    if (window_params->sim_running) {
        ////////////////////////////////////////////////////////////////
        // calculate forces between all body pairs
        ////////////////////////////////////////////////////////////////
        if (global_bodies->count > 0) {
            const double delta_t= window_params->time_step;

            // Step 1: update all positions
            // x(t + Δt) = x(t) + v(t) * Δt + 0.5 * a(t) * Δt^2
            for (int i = 0; i < global_bodies->count; i++) {
                body_t* body = &global_bodies->bodies[i];
                const vec3 pos_delta = vec3_add(
                    vec3_scale(body->vel, delta_t),
                    vec3_scale(body->acc, 0.5 *delta_t* delta_t)
                );
                body->pos = vec3_add(body->pos, pos_delta);
                body->acc_prev = body->acc;
            }

            // Step 2: calculate all forces using new positions
            for (int i = 0; i < global_bodies->count; i++) {
                global_bodies->bodies[i].force = vec3_zero();
            }
            for (int i = 0; i < global_bodies->count; i++) {
                for (int j = i + 1; j < global_bodies->count; j++) {
                    body_t* force_receiving_body = &global_bodies->bodies[i];
                    body_t* force_applying_body = &global_bodies->bodies[j];

                    // third law pair: the force is applied to both bodies
                    const vec3 grav_force =  body_calculateGravForce(sim, i, j);
                    force_receiving_body->force = vec3_add(force_receiving_body->force, grav_force);
                    force_applying_body->force = vec3_sub(force_applying_body->force, grav_force);
                }
            }

            // Step 3: update accelerations and velocities
            // v(t + Δt) = v(t) + 0.5 * (a(t) + a(t + Δt)) * Δt
            for (int i = 0; i < global_bodies->count; i++) {
                body_t* body = &global_bodies->bodies[i];
                body->acc = vec3_scale(body->force, 1.0 / body->mass);
                const vec3 avg_acc = vec3_add(body->acc_prev, body->acc);
                body->vel = vec3_add(body->vel, vec3_scale(avg_acc, 0.5 * delta_t));
                body->vel_mag = vec3_mag(body->vel);

                body_updateRotation(body, delta_t);
            }
        }

        ////////////////////////////////////////////////////////////////
        // calculate forces between spacecraft and bodies
        ////////////////////////////////////////////////////////////////
        if (global_spacecraft->count > 0 && global_bodies->count > 0) {
            const double delta_t = window_params->time_step;

            // Step 1: update all spacecraft positions
            // x(t + Δt) = x(t) + v(t) * Δt + 0.5 * a(t) * Δt^2
            for (int i = 0; i < global_spacecraft->count; i++) {
                spacecraft_t* craft = &global_spacecraft->spacecraft[i];
                const vec3 pos_delta = vec3_add(
                    vec3_scale(craft->vel, delta_t),
                    vec3_scale(craft->acc, 0.5 *delta_t* delta_t)
                );
                craft->pos = vec3_add(craft->pos, pos_delta);
                craft->acc_prev = craft->acc;
            }

            // Step 2: calculate all forces using new positions
            for (int i = 0; i < global_spacecraft->count; i++) {
                spacecraft_t* craft = &global_spacecraft->spacecraft[i];
                craft->grav_force = vec3_zero();

                // check if burn should be active
                craft_checkBurnSchedule(craft, global_bodies, window_params->sim_time);

                // calculate gravitational forces from all bodies
                for (int j = 0; j < global_bodies->count; j++) {
                    craft->current_total_mass = craft->fuel_mass + craft->dry_mass;
                    craft->grav_force = vec3_add(craft->grav_force, craft_calculateGravForce(sim, i, j));
                }
                // add J2 perturbation force to craft
                // apply J2 perturbation force for the closest body
                craft->grav_force = vec3_add(craft->grav_force, craft_calculateJ2Force(sim, i, craft->orbital_elements.SOI_planet_id));

                // apply thrust and consume fuel
                craft_applyThrust(craft);
                craft_consumeFuel(craft, delta_t);
            }

            // Step 3: update accelerations and velocities
            // v(t + Δt) = v(t) + 0.5 * (a(t) + a(t + Δt)) * Δt
            for (int i = 0; i < global_spacecraft->count; i++) {
                spacecraft_t* craft = &global_spacecraft->spacecraft[i];
                craft->acc = vec3_scale(craft->grav_force, 1.0 / craft->current_total_mass);
                const vec3 avg_acc = vec3_add(craft->acc_prev, craft->acc);
                craft->vel = vec3_add(craft->vel, vec3_scale(avg_acc, 0.5 * delta_t));
                craft->vel_mag = vec3_mag(craft->vel);

                // calculate orbital elements relative to the SOI body (or closest body)
                if (craft->orbital_elements.SOI_planet_id >= 0 && craft->orbital_elements.SOI_planet_id < global_bodies->count) {
                    // while the orbital elements for bodies are calculated outside the sim thread, the spacecraft
                    // elements are so the craft is "aware" of its actual position asap this is potentially useful
                    // if the craft were to burn automatically to try and achieve a certain orbital parameter
                    calculateOrbitalElements(&craft->orbital_elements, &craft->pos, &craft->vel, &global_bodies->bodies[craft->orbital_elements.SOI_planet_id]);
                }
            }
        }

        // increment simulation time
        if (global_bodies->count > 0) {
            window_params->sim_time += window_params->time_step;
        }
    }
}