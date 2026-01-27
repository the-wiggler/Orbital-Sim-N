#include "simulation.h"
#include "../globals.h"
#include "../sim/bodies.h"
#include "../sim/spacecraft.h"
#include "../math/matrix.h"
#include <math.h>
#include <stdlib.h>

// calculate total system energy for all bodies
double calculateTotalSystemEnergy(const sim_properties_t* sim) {
    const body_properties_t* gb = &sim->gb;
    const spacecraft_properties_t* sc = &sim->gs;

    double total_kinetic = 0.0;
    double total_potential = 0.0;

    // calculate kinetic energy for all bodies
    for (int i = 0; i < gb->count; i++) {
        const body_t* body = &gb->bodies[i];
        total_kinetic += 0.5 * body->mass * body->vel_mag * body->vel_mag;
    }

    // calculate kinetic energy for all spacecraft
    for (int i = 0; i < sc->count; i++) {
        const spacecraft_t* craft = &sc->spacecraft[i];
        total_kinetic += 0.5 * craft->current_total_mass * craft->vel_mag * craft->vel_mag;
    }

    // calculate potential energy between all body pairs
    for (int i = 0; i < gb->count; i++) {
        for (int j = i + 1; j < gb->count; j++) {
            const body_t* bi = &gb->bodies[i];
            const body_t* bj = &gb->bodies[j];
            const vec3 delta = vec3_sub(bj->pos, bi->pos);
            const double r = vec3_mag(delta);
            if (r > 0) {
                total_potential += -(G * bi->mass * bj->mass) / r;
            }
        }
    }

    // calculate potential energy between spacecraft and bodies
    for (int i = 0; i < sc->count; i++) {
        for (int j = 0; j < gb->count; j++) {
            const spacecraft_t* craft = &sc->spacecraft[i];
            const body_t* body = &gb->bodies[j];
            const vec3 delta = vec3_sub(body->pos, craft->pos);
            const double r = vec3_mag(delta);
            if (r > 0) {
                total_potential += -(G * craft->current_total_mass * body->mass) / r;
            }
        }
    }

    return total_kinetic + total_potential;
}

// calculates orbital elements (this probably needs to be optimized somehow at some point because this seems very resource heavy)
void calculateOrbitalElements(orbital_elements_t* target_orbital_elements, vec3 target_pos, vec3 target_vel, const body_t* body) {
    orbital_elements_t* oe = target_orbital_elements;
    // first, the initial properties of the craft relative to the target planet should be calculated
    const vec3 c_pos     = vec3_sub(target_pos, body->pos); // position vector
    const vec3 c_vel     = vec3_sub(target_vel, body->vel); // velocity vector
    const double c_r     = vec3_mag(c_pos); // distance
    const double c_speed = vec3_mag(c_vel);
    const double mu      = G * body->mass; // gravitational parameter
    const vec3 c_h       = vec3_cross(c_pos, c_vel); // specific angular momentum
    const vec3 k         = { 0, 0, 1 };
    const vec3 c_n       = vec3_cross(k, c_h); // ascending node vector

    const vec3 term1     = vec3_scalar_div(vec3_cross(c_vel, c_h), mu); // first term of e_vec
    const vec3 term2     = vec3_scalar_div(c_pos, c_r); // second term of e_vec
    const vec3 e_vec     = vec3_sub(term1, term2); // eccentricity vector
    oe->specific_E    = ((c_speed * c_speed) / 2) - (mu / c_r); // specific orbital energy

    // orbital elements
    oe->semi_major_axis = -1.0 * (mu / (2 * oe->specific_E));
    oe->eccentricity = vec3_mag(e_vec);

    const double h_mag = vec3_mag(c_h);
    oe->inclination = acos(c_h.z / h_mag); // the angle between the orbital and equatorial planes

    // longitude of ascending node -- the angle from the vernal equinox vector to the ascending node on the equatorial plane
    const double n_mag = vec3_mag(c_n);
    if (n_mag > 1e-10) {
        oe->ascending_node = atan2(c_n.y, c_n.x);
        if (oe->ascending_node < 0) {
            oe->ascending_node += 2 * PI;
        }
    } else {
        oe->ascending_node = 0.0; // undefined for equatorial orbits (probably unlikely to happen perfectly)
    }

    // argument of periapsis -- the angle measured between the ascending node and the perigee
    if (oe->eccentricity > 1e-10 && n_mag > 1e-10) {
        const double cos_omega = vec3_dot(c_n, e_vec) / (n_mag * oe->eccentricity);
        oe->arg_periapsis = acos(fmax(-1.0, fmin(1.0, cos_omega)));
        if (e_vec.z < 0) {
            oe->arg_periapsis = 2 * PI - oe->arg_periapsis;
        }
    } else if (oe->eccentricity > 1e-10) {
        // equatorial orbit, use longitude of periapsis
        oe->arg_periapsis = atan2(e_vec.y, e_vec.x);
        if (oe->arg_periapsis < 0) {
            oe->arg_periapsis += 2 * PI;
        }
    } else {
        oe->arg_periapsis = 0.0; // undefined for circular orbits
    }

    // true anomaly -- the angle between perigee and satellite in the orbital plane at a specific time
    if (oe->eccentricity > 1e-10) {
        const double cos_nu = vec3_dot(e_vec, c_pos) / (oe->eccentricity * c_r);
        oe->true_anomaly = acos(fmax(-1.0, fmin(1.0, cos_nu)));
        if (vec3_dot(c_pos, c_vel) < 0) {
            oe->true_anomaly = 2 * PI - oe->true_anomaly;
        }
    } else {
        // circular orbit, use argument of latitude
        if (n_mag > 1e-10) {
            const double cos_u = vec3_dot(c_n, c_pos) / (n_mag * c_r);
            oe->true_anomaly = acos(fmax(-1.0, fmin(1.0, cos_u)));
            if (c_pos.z < 0) {
                oe->true_anomaly = 2 * PI - oe->true_anomaly;
            }
        } else {
            // equatorial and circular, use true longitude
            oe->true_anomaly = atan2(c_pos.y, c_pos.x);
            if (oe->true_anomaly < 0) {
                oe->true_anomaly += 2 * PI;
            }
        }
    }
}

// reset the simulation by removing all bodies from the system
void resetSim(sim_properties_t* sim) {
    body_properties_t* gb = &sim->gb;
    spacecraft_properties_t* sc = &sim->gs;
    window_params_t* wp = &sim->wp;

    // reset simulation time
    wp->sim_time = 0;
    wp->reset_sim = false;

    // reset all bodies
    gb->count = 0;

    // reset all spacecraft
    sc->count = 0;
}


// probably the most important function in the entire program lol. This is responsible for updating the position of literally
// every single planet each time step. It uses a method called Velocity Verlet Integration.
// https://en.wikipedia.org/wiki/Verlet_integration#Velocity_Verlet
void runCalculations(sim_properties_t* sim) {
    body_properties_t* gb = &sim->gb;
    spacecraft_properties_t* sc = &sim->gs;
    window_params_t* wp = &sim->wp;

    if (wp->sim_running) {
        ////////////////////////////////////////////////////////////////
        // calculate forces between all body pairs
        ////////////////////////////////////////////////////////////////
        if (gb->count > 0) {
            const double dt = wp->time_step;

            // Step 1: update all positions
            // x(t + Δt) = x(t) + v(t) * Δt + 0.5 * a(t) * Δt^2
            for (int i = 0; i < gb->count; i++) {
                body_t* body = &gb->bodies[i];
                const vec3 pos_delta = vec3_add(
                    vec3_scale(body->vel, dt),
                    vec3_scale(body->acc, 0.5 * dt * dt)
                );
                body->pos = vec3_add(body->pos, pos_delta);
                body->acc_prev = body->acc;
            }

            // Step 2: calculate all forces using new positions
            for (int i = 0; i < gb->count; i++) {
                gb->bodies[i].force = vec3_zero();
                gb->bodies[i].oe.closest_r_squared = INFINITY;
            }
            for (int i = 0; i < gb->count; i++) {
                for (int j = i + 1; j < gb->count; j++) {
                    body_calculateGravForce(sim, i, j);
                }
            }

            // check if bodies have exited their current SOI
            for (int i = 0; i < gb->count; i++) {
                body_t* body = &gb->bodies[i];
                if (body->oe.SOI_planet_id > 0 && body->oe.SOI_planet_id < gb->count) {
                    const body_t* soi_body = &gb->bodies[body->oe.SOI_planet_id];
                    const vec3 delta = vec3_sub(soi_body->pos, body->pos);
                    const double dist = vec3_mag(delta);
                    if (dist > soi_body->SOI_radius) {
                        // exited SOI, fall back to the closest body
                        body->oe.SOI_planet_id = body->oe.closest_planet_id;
                    }
                }
            }

            // Step 3: update accelerations and velocities
            // v(t + Δt) = v(t) + 0.5 * (a(t) + a(t + Δt)) * Δt
            for (int i = 0; i < gb->count; i++) {
                body_t* body = &gb->bodies[i];
                body->acc = vec3_scale(body->force, 1.0 / body->mass);
                const vec3 avg_acc = vec3_add(body->acc_prev, body->acc);
                body->vel = vec3_add(body->vel, vec3_scale(avg_acc, 0.5 * dt));
                body->vel_mag = vec3_mag(body->vel);

                body_calculateKineticEnergy(body);
                body_updateRotation(body, dt);
            }

            // calculate orbital elements for bodies (skip central body at index 0)
            for (int i = 1; i < gb->count; i++) {
                body_t* body = &gb->bodies[i];
                if (body->oe.SOI_planet_id != i && body->oe.SOI_planet_id >= 0 && body->oe.SOI_planet_id < gb->count) {
                    calculateOrbitalElements(&body->oe, body->pos, body->vel, &gb->bodies[body->oe.SOI_planet_id]);
                }
            }
        }

        ////////////////////////////////////////////////////////////////
        // calculate forces between spacecraft and bodies :)
        ////////////////////////////////////////////////////////////////
        if (sc->count > 0 && gb->count > 0) {
            const double dt = wp->time_step;

            // Step 1: update all spacecraft positions
            // x(t + Δt) = x(t) + v(t) * Δt + 0.5 * a(t) * Δt^2
            for (int i = 0; i < sc->count; i++) {
                spacecraft_t* craft = &sc->spacecraft[i];
                const vec3 pos_delta = vec3_add(
                    vec3_scale(craft->vel, dt),
                    vec3_scale(craft->acc, 0.5 * dt * dt)
                );
                craft->pos = vec3_add(craft->pos, pos_delta);
                craft->acc_prev = craft->acc;
            }

            // Step 2: calculate all forces using new positions
            for (int i = 0; i < sc->count; i++) {
                spacecraft_t* craft = &sc->spacecraft[i];
                craft->grav_force = vec3_zero();
                craft->oe.closest_r_squared = INFINITY;

                // check if burn should be active
                craft_checkBurnSchedule(craft, gb, wp->sim_time);

                // calculate gravitational forces from all bodies
                for (int j = 0; j < gb->count; j++) {
                    craft_calculateGravForce(sim, i, j);
                }

                // check if spacecraft has exited its current SOI
                if (craft->oe.SOI_planet_id > 0 && craft->oe.SOI_planet_id < gb->count) {
                    const body_t* soi_body = &gb->bodies[craft->oe.SOI_planet_id];
                    const vec3 delta = vec3_sub(soi_body->pos, craft->pos);
                    const double dist = vec3_mag(delta);
                    if (dist > soi_body->SOI_radius) {
                        // exited SOI, fall back to the closest body
                        craft->oe.SOI_planet_id = craft->oe.closest_planet_id;
                    }
                }

                // apply thrust and consume fuel
                craft_applyThrust(craft);
                craft_consumeFuel(craft, dt);
            }

            // Step 3: update accelerations and velocities
            // v(t + Δt) = v(t) + 0.5 * (a(t) + a(t + Δt)) * Δt
            for (int i = 0; i < sc->count; i++) {
                spacecraft_t* craft = &sc->spacecraft[i];
                craft->acc = vec3_scale(craft->grav_force, 1.0 / craft->current_total_mass);
                const vec3 avg_acc = vec3_add(craft->acc_prev, craft->acc);
                craft->vel = vec3_add(craft->vel, vec3_scale(avg_acc, 0.5 * dt));
                craft->vel_mag = vec3_mag(craft->vel);

                // calculate orbital elements relative to the SOI body (or closest body)
                if (craft->oe.SOI_planet_id >= 0 && craft->oe.SOI_planet_id < gb->count) {
                    calculateOrbitalElements(&craft->oe, craft->pos, craft->vel, &gb->bodies[craft->oe.SOI_planet_id]);
                }
            }
        }

        // increment simulation time
        if (gb->count > 0) {
            wp->sim_time += wp->time_step;
        }
    }
}

// cleanup for main
void cleanup(const sim_properties_t* sim) {
    // no memory cleanup needed with static allocation
    (void)sim;  // suppress unused parameter warning
}
