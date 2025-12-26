#include "simulation.h"
#include "../globals.h"
#include "../sim/bodies.h"
#include "../sim/spacecraft.h"
#include <math.h>
#include <stdlib.h>

// calculate total system energy for all bodies
// this avoids double-counting by only calculating each pair interaction once
double calculateTotalSystemEnergy(const sim_properties_t* sim) {
    const body_properties_t* gb = &sim->gb;
    const spacecraft_properties_t* sc = &sim->gs;

    double total_kinetic = 0.0;
    double total_potential = 0.0;

    // calculate kinetic energy for all bodies
    for (int i = 0; i < gb->count; i++) {
        total_kinetic += 0.5 * gb->mass[i] * gb->vel[i] * gb->vel[i];
    }

    // calculate kinetic energy for all spacecraft
    for (int i = 0; i < sc->count; i++) {
        total_kinetic += 0.5 * sc->current_total_mass[i] * sc->vel[i] * sc->vel[i];
    }

    // calculate potential energy between all body pairs
    for (int i = 0; i < gb->count; i++) {
        for (int j = i + 1; j < gb->count; j++) {
            double delta_x = gb->pos_x[j] - gb->pos_x[i];
            double delta_y = gb->pos_y[j] - gb->pos_y[i];
            double r = sqrt(delta_x * delta_x + delta_y * delta_y);
            if (r > 0) {  // prevent divide by zero
                total_potential += -(G * gb->mass[i] * gb->mass[j]) / r;
            }
        }
    }

    // calculate potential energy between spacecraft and bodies
    for (int i = 0; i < sc->count; i++) {
        for (int j = 0; j < gb->count; j++) {
            double delta_x = gb->pos_x[j] - sc->pos_x[i];
            double delta_y = gb->pos_y[j] - sc->pos_y[i];
            double r = sqrt(delta_x * delta_x + delta_y * delta_y);
            if (r > 0) {  // prevent divide by zero
                total_potential += -(G * sc->current_total_mass[i] * gb->mass[j]) / r;
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

    // reset simulation time to 0
    wp->sim_time = 0;
    // change sim reset flag back to false
    wp->reset_sim = false;

    // free all bodies from memory
    if (gb != NULL) {
        for (int i = 0; i < gb->count; i++) {
            free(gb->names[i]);
        }

        // free cache arrays for each body
        if (gb->cached_body_coords_x != NULL) {
            free(gb->cached_body_coords_x);
        }
        if (gb->cached_body_coords_y != NULL) {
            free(gb->cached_body_coords_y);
        }

        free(gb->names);
        free(gb->mass);
        free(gb->radius);
        free(gb->pixel_radius);
        free(gb->pos_x);
        free(gb->pos_y);
        free(gb->pos_x);
        free(gb->pos_y);
        free(gb->vel_x);
        free(gb->vel_y);
        free(gb->vel);
        free(gb->acc_x);
        free(gb->acc_y);
        free(gb->acc_x_prev);
        free(gb->acc_y_prev);
        free(gb->force_x);
        free(gb->force_y);
        free(gb->kinetic_energy);

        gb->count = 0;
        gb->names = NULL;
        gb->mass = NULL;
        gb->radius = NULL;
        gb->pixel_radius = NULL;
        gb->ndc_x = NULL;
        gb->ndc_y = NULL;
        gb->pos_x = NULL;
        gb->pos_y = NULL;
        gb->vel_x = NULL;
        gb->vel_y = NULL;
        gb->vel = NULL;
        gb->acc_x = NULL;
        gb->acc_y = NULL;
        gb->acc_x_prev = NULL;
        gb->acc_y_prev = NULL;
        gb->force_x = NULL;
        gb->force_y = NULL;
        gb->kinetic_energy = NULL;
        gb->cached_body_coords_x = NULL;
        gb->cached_body_coords_y = NULL;
        gb->cache_counter = 0;
        gb->cache_valid_count = 0;
    }

    // free all craft from memory
    if (sc != NULL) {
        for (int i = 0; i < sc->count; i++) {
            free(sc->names[i]);
            free(sc->burn_properties[i]);
        }
        free(sc->names);
        free(sc->current_total_mass);
        free(sc->dry_mass);
        free(sc->fuel_mass);
        free(sc->pos_x);
        free(sc->pos_y);
        free(sc->attitude);
        free(sc->vel_x);
        free(sc->vel_y);
        free(sc->vel);
        free(sc->rotational_v);
        free(sc->momentum);
        free(sc->acc_x);
        free(sc->acc_y);
        free(sc->acc_x_prev);
        free(sc->acc_y_prev);
        free(sc->rotational_a);
        free(sc->moment_of_inertia);
        free(sc->grav_force_x);
        free(sc->grav_force_y);
        free(sc->torque);
        free(sc->thrust);
        free(sc->mass_flow_rate);
        free(sc->specific_impulse);
        free(sc->throttle);
        free(sc->nozzle_gimbal_range);
        free(sc->nozzle_velocity);
        free(sc->engine_on);
        free(sc->num_burns);
        free(sc->burn_properties);

        sc->count = 0;
        sc->names = NULL;
        sc->current_total_mass = NULL;
        sc->dry_mass = NULL;
        sc->fuel_mass = NULL;
        sc->pos_x = NULL;
        sc->pos_y = NULL;
        sc->attitude = NULL;
        sc->vel_x = NULL;
        sc->vel_y = NULL;
        sc->vel = NULL;
        sc->rotational_v = NULL;
        sc->momentum = NULL;
        sc->acc_x = NULL;
        sc->acc_y = NULL;
        sc->acc_x_prev = NULL;
        sc->acc_y_prev = NULL;
        sc->rotational_a = NULL;
        sc->moment_of_inertia = NULL;
        sc->grav_force_x = NULL;
        sc->grav_force_y = NULL;
        sc->torque = NULL;
        sc->thrust = NULL;
        sc->mass_flow_rate = NULL;
        sc->specific_impulse = NULL;
        sc->throttle = NULL;
        sc->nozzle_gimbal_range = NULL;
        sc->nozzle_velocity = NULL;
        sc->engine_on = NULL;
        sc->num_burns = NULL;
        sc->burn_properties = NULL;
    }
}

void runCalculations(sim_properties_t* sim) {
    body_properties_t* gb = &sim->gb;
    spacecraft_properties_t* sc = &sim->gs;
    window_params_t* wp = &sim->wp;

    if (wp->sim_running) {
        ////////////////////////////////////////////////////////////////
        // calculate forces between all body pairs
        ////////////////////////////////////////////////////////////////
        if (gb != NULL && gb->count > 0) {
            // initialize forces to zero to re-calculate them
            for (int i = 0; i < gb->count; i++) {
                gb->force_x[i] = 0;
                gb->force_y[i] = 0;
            }

            // loop through every body and add the resultant force to the subject body force vector
            for (int i = 0; i < gb->count; i++) {
                for (int j = i + 1; j < gb->count; j++) {
                    body_calculateGravForce(sim, i, j);
                }
            }

            // calculate kinetic energy for each body
            for (int i = 0; i < gb->count; i++) {
                body_calculateKineticEnergy(gb, i);
            }

            // update the motion for each body and draw
            for (int i = 0; i < gb->count; i++) {
                // updates the kinematic properties of each body (velocity, acceleration, position, etc.)
                body_updateMotion(gb, i, wp->time_step);
            }
        }

        ////////////////////////////////////////////////////////////////
        // calculate forces between spacecraft and bodies
        ////////////////////////////////////////////////////////////////
        if (sc != NULL && sc->count > 0 && gb != NULL && gb->count > 0) {
            for (int i = 0; i < sc->count; i++) {
                sc->grav_force_x[i] = 0;
                sc->grav_force_y[i] = 0;

                // check if burn should be active based on simulation time
                craft_checkBurnSchedule(sc, i, gb, wp->sim_time);

                // loop through all bodies and calculate gravitational forces on spacecraft
                for (int j = 0; j < gb->count; j++) {
                    craft_calculateGravForce(sim, i, j);
                }

                // apply thrust first, then consume fuel
                craft_applyThrust(sc, i);
                craft_consumeFuel(sc, i, wp->time_step);
            }

            // update the motion for each craft and draw
            for (int i = 0; i < sc->count; i++) {
                // updates the kinematic properties of each body (velocity, acceleration, position, etc.)
                craft_updateMotion(sc, i, wp->time_step);
            }
        }

        // increment the time based on the time step
        if (gb != NULL && gb->count > 0) {
            wp->sim_time += wp->time_step;
        }
    }
}

// cleanup for main
void cleanup(sim_properties_t* sim) {
    body_properties_t* gb = &sim->gb;
    spacecraft_properties_t* sc = &sim->gs;

    // free all bodies
    for (int i = 0; i < gb->count; i++) {
        free(gb->names[i]);
    }

    // free cache arrays for each body
    if (gb->cached_body_coords_x != NULL) {
        for (int i = 0; i < gb->count; i++) {
            free(gb->cached_body_coords_x[i]);
        }
        free(gb->cached_body_coords_x);
    }
    if (gb->cached_body_coords_y != NULL) {
        for (int i = 0; i < gb->count; i++) {
            free(gb->cached_body_coords_y[i]);
        }
        free(gb->cached_body_coords_y);
    }

    free(gb->names);
    free(gb->mass);
    free(gb->radius);
    free(gb->pixel_radius);
    free(gb->pos_x);
    free(gb->pos_y);
    free(gb->vel_x);
    free(gb->vel_y);
    free(gb->vel);
    free(gb->acc_x);
    free(gb->acc_y);
    free(gb->acc_x_prev);
    free(gb->acc_y_prev);
    free(gb->force_x);
    free(gb->force_y);
    free(gb->kinetic_energy);

    // free all spacecraft
    for (int i = 0; i < sc->count; i++) {
        free(sc->names[i]);
        free(sc->burn_properties[i]);
    }
    free(sc->names);
    free(sc->current_total_mass);
    free(sc->dry_mass);
    free(sc->fuel_mass);
    free(sc->pos_x);
    free(sc->pos_y);
    free(sc->attitude);
    free(sc->vel_x);
    free(sc->vel_y);
    free(sc->vel);
    free(sc->rotational_v);
    free(sc->momentum);
    free(sc->acc_x);
    free(sc->acc_y);
    free(sc->acc_x_prev);
    free(sc->acc_y_prev);
    free(sc->rotational_a);
    free(sc->moment_of_inertia);
    free(sc->grav_force_x);
    free(sc->grav_force_y);
    free(sc->torque);
    free(sc->thrust);
    free(sc->mass_flow_rate);
    free(sc->specific_impulse);
    free(sc->throttle);
    free(sc->nozzle_gimbal_range);
    free(sc->nozzle_velocity);
    free(sc->engine_on);
    free(sc->num_burns);
    free(sc->burn_properties);
}
