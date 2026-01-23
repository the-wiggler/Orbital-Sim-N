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

// reset the simulation by removing all bodies from the system
void resetSim(sim_properties_t* sim) {
    body_properties_t* gb = &sim->gb;
    spacecraft_properties_t* sc = &sim->gs;
    window_params_t* wp = &sim->wp;

    // reset simulation time
    wp->sim_time = 0;
    wp->reset_sim = false;

    // free all bodies
    if (gb->bodies != NULL) {
        for (int i = 0; i < gb->count; i++) {
            free(gb->bodies[i].name);
        }
        free(gb->bodies);
        gb->bodies = NULL;
        gb->count = 0;
        gb->capacity = 0;
    }

    // free all spacecraft
    if (sc->spacecraft != NULL) {
        for (int i = 0; i < sc->count; i++) {
            free(sc->spacecraft[i].name);
            free(sc->spacecraft[i].burn_properties);
        }
        free(sc->spacecraft);
        sc->spacecraft = NULL;
        sc->count = 0;
        sc->capacity = 0;
    }
}


// probably the most important function in the entire program lol. This is responsible for updating the position of literally
// every single planet each time step. It uses a method called Velocity Verlet Integration.
// https://en.wikipedia.org/wiki/Verlet_integration#Velocity_Verlet
void runCalculations(sim_properties_t* sim) {
    const body_properties_t* gb = &sim->gb;
    const spacecraft_properties_t* sc = &sim->gs;
    window_params_t* wp = &sim->wp;

    if (wp->sim_running) {
        ////////////////////////////////////////////////////////////////
        // calculate forces between all body pairs
        ////////////////////////////////////////////////////////////////
        if (gb->bodies != NULL && gb->count > 0) {
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
            }
            for (int i = 0; i < gb->count; i++) {
                for (int j = i + 1; j < gb->count; j++) {
                    body_calculateGravForce(sim, i, j);
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
        }

        ////////////////////////////////////////////////////////////////
        // calculate forces between spacecraft and bodies :)
        ////////////////////////////////////////////////////////////////
        if (sc->spacecraft != NULL && sc->count > 0 && gb->bodies != NULL && gb->count > 0) {
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
                craft->closest_r_squared = INFINITY;

                // check if burn should be active
                craft_checkBurnSchedule(craft, gb, wp->sim_time);

                // calculate gravitational forces from all bodies
                for (int j = 0; j < gb->count; j++) {
                    craft_calculateGravForce(sim, i, j);
                }

                // check if spacecraft has exited its current SOI
                if (craft->SOI_planet_id > 0 && craft->SOI_planet_id < gb->count) {
                    const body_t* soi_body = &gb->bodies[craft->SOI_planet_id];
                    const vec3 delta = vec3_sub(soi_body->pos, craft->pos);
                    const double dist = vec3_mag(delta);
                    if (dist > soi_body->SOI_radius) {
                        // exited SOI, fall back to the closest body
                        craft->SOI_planet_id = craft->closest_planet_id;
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
                if (craft->SOI_planet_id >= 0 && craft->SOI_planet_id < gb->count) {
                    craft_calculateOrbitalElements(craft, &gb->bodies[craft->SOI_planet_id]);
                }
            }
        }

        // increment simulation time
        if (gb->bodies != NULL && gb->count > 0) {
            wp->sim_time += wp->time_step;
        }
    }
}

// cleanup for main
void cleanup(const sim_properties_t* sim) {
    const body_properties_t* gb = &sim->gb;
    const spacecraft_properties_t* sc = &sim->gs;

    // free all bodies
    if (gb->bodies != NULL) {
        for (int i = 0; i < gb->count; i++) {
            free(gb->bodies[i].name);
        }
        free(gb->bodies);
    }

    // free all spacecraft
    if (sc->spacecraft != NULL) {
        for (int i = 0; i < sc->count; i++) {
            free(sc->spacecraft[i].name);
            free(sc->spacecraft[i].burn_properties);
        }
        free(sc->spacecraft);
    }
}
