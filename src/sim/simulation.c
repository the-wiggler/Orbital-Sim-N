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

void runCalculations(sim_properties_t* sim) {
    const body_properties_t* gb = &sim->gb;
    const spacecraft_properties_t* sc = &sim->gs;
    window_params_t* wp = &sim->wp;

    if (wp->sim_running) {
        ////////////////////////////////////////////////////////////////
        // calculate forces between all body pairs
        ////////////////////////////////////////////////////////////////
        if (gb->bodies != NULL && gb->count > 0) {
            // reset forces to zero
            for (int i = 0; i < gb->count; i++) {
                gb->bodies[i].force = vec3_zero();
            }

            // calculate gravitational forces between all body pairs
            for (int i = 0; i < gb->count; i++) {
                for (int j = i + 1; j < gb->count; j++) {
                    body_calculateGravForce(sim, i, j);
                }
            }

            // calculate kinetic energy and update motion for each body
            for (int i = 0; i < gb->count; i++) {
                body_t* body = &gb->bodies[i];
                body_calculateKineticEnergy(body);
                body_updateMotion(body, wp->time_step);
                body_updateRotation(body, wp->time_step);
            }
        }

        ////////////////////////////////////////////////////////////////
        // calculate forces between spacecraft and bodies
        ////////////////////////////////////////////////////////////////
        if (sc->spacecraft != NULL && sc->count > 0 && gb->bodies != NULL && gb->count > 0) {
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
                craft_consumeFuel(craft, wp->time_step);
            }

            // update motion and orbital elements for each craft
            for (int i = 0; i < sc->count; i++) {
                spacecraft_t* craft = &sc->spacecraft[i];
                craft_updateMotion(craft, wp->time_step);

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
