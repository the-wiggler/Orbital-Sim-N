#include "spacecraft.h"
#include "../globals.h"
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include "../math/matrix.h"

void displayError(const char* title, const char* message);

// calculates orbital elements (this probably needs to be optimized somehow at some point because this seems very resource heavy)
void craft_calculateOrbitalElements(spacecraft_t* craft, const body_t* body) {
    // first, the initial properties of the craft relative to the target planet should be calculated
    const vec3 c_pos     = vec3_sub(craft->pos, body->pos); // position vector
    const vec3 c_vel     = vec3_sub(craft->vel, body->vel); // velocity vector
    const double c_r     = vec3_mag(c_pos); // distance
    const double c_speed = vec3_mag(c_vel);
    const double mu      = G * body->mass; // gravitational parameter
    const vec3 c_h       = vec3_cross(c_pos, c_vel); // specific angular momentum
    const vec3 k         = { 0, 0, 1 };
    const vec3 c_n       = vec3_cross(k, c_h); // ascending node vector

    const vec3 term1     = vec3_scalar_div(vec3_cross(c_vel, c_h), mu); // first term of e_vec
    const vec3 term2     = vec3_scalar_div(c_pos, c_r); // second term of e_vec
    const vec3 e_vec     = vec3_sub(term1, term2); // eccentricity vector
    craft->specific_E    = ((c_speed * c_speed) / 2) - (mu / c_r); // specific orbital energy

    // orbital elements
    craft->semi_major_axis = -1.0 * (mu / (2 * craft->specific_E));
    craft->eccentricity = vec3_mag(e_vec);

    const double h_mag = vec3_mag(c_h);
    craft->inclination = acos(c_h.z / h_mag); // the angle between the orbital and equatorial planes

    // longitude of ascending node -- the angle from the vernal equinox vector to the ascending node on the equatorial plane
    const double n_mag = vec3_mag(c_n);
    if (n_mag > 1e-10) {
        craft->ascending_node = atan2(c_n.y, c_n.x);
        if (craft->ascending_node < 0) {
            craft->ascending_node += 2 * PI;
        }
    } else {
        craft->ascending_node = 0.0; // undefined for equatorial orbits (probably unlikely to happen perfectly)
    }

    // argument of periapsis -- the angle measured between the ascending node and the perigee
    if (craft->eccentricity > 1e-10 && n_mag > 1e-10) {
        const double cos_omega = vec3_dot(c_n, e_vec) / (n_mag * craft->eccentricity);
        craft->arg_periapsis = acos(fmax(-1.0, fmin(1.0, cos_omega)));
        if (e_vec.z < 0) {
            craft->arg_periapsis = 2 * PI - craft->arg_periapsis;
        }
    } else if (craft->eccentricity > 1e-10) {
        // equatorial orbit, use longitude of periapsis
        craft->arg_periapsis = atan2(e_vec.y, e_vec.x);
        if (craft->arg_periapsis < 0) {
            craft->arg_periapsis += 2 * PI;
        }
    } else {
        craft->arg_periapsis = 0.0; // undefined for circular orbits
    }

    // true anomaly -- the angle between perigee and satellite in the orbital plane at a specific time
    if (craft->eccentricity > 1e-10) {
        const double cos_nu = vec3_dot(e_vec, c_pos) / (craft->eccentricity * c_r);
        craft->true_anomaly = acos(fmax(-1.0, fmin(1.0, cos_nu)));
        if (vec3_dot(c_pos, c_vel) < 0) {
            craft->true_anomaly = 2 * PI - craft->true_anomaly;
        }
    } else {
        // circular orbit, use argument of latitude
        if (n_mag > 1e-10) {
            const double cos_u = vec3_dot(c_n, c_pos) / (n_mag * c_r);
            craft->true_anomaly = acos(fmax(-1.0, fmin(1.0, cos_u)));
            if (c_pos.z < 0) {
                craft->true_anomaly = 2 * PI - craft->true_anomaly;
            }
        } else {
            // equatorial and circular, use true longitude
            craft->true_anomaly = atan2(c_pos.y, c_pos.x);
            if (craft->true_anomaly < 0) {
                craft->true_anomaly += 2 * PI;
            }
        }
    }
}

// check and activate burns
void craft_checkBurnSchedule(spacecraft_t* craft, const body_properties_t* gb, const double sim_time) {
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
            const body_t* target = &gb->bodies[target_id];

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
void craft_calculateGravForce(sim_properties_t* sim, const int craft_idx, const int body_idx) {
    spacecraft_t* craft = &sim->gs.spacecraft[craft_idx];
    const body_t* body = &sim->gb.bodies[body_idx];

    // calculate the distance between the spacecraft and the body
    const vec3 delta_pos = vec3_sub(body->pos, craft->pos);
    const double r_squared = vec3_mag_sq(delta_pos);
    const double r = sqrt(r_squared);

    // planet collision logic
    if (r < body->radius) {
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
    const double r_cubed = r_squared * r;
    const double force_factor = (G * craft->current_total_mass * body->mass) / r_cubed;

    // apply the force to the craft
    const vec3 force = vec3_scale(delta_pos, force_factor);
    craft->grav_force = vec3_add(craft->grav_force, force);

    // checks for new closest planet and checks if its in the SOI
    if (r_squared < craft->closest_r_squared) {
        craft->closest_r_squared = r_squared;
        craft->closest_planet_id = body_idx;
        if (r <= body->SOI_radius) {
            craft->SOI_planet_id = body_idx;
        }
    }
}

// updates the ID and distance of the closest planet
// (this should be used when initially spawning a craft because the grav calculations do this exact calculation by default)
// this is probably executed when the JSON is loaded
void craft_findClosestPlanet(spacecraft_t* craft, const body_properties_t* gb) {
    double closest_r_squared = INFINITY;
    int closest_planet_id = 0;
    for (int i = 0; i < gb->count; i++) {
        // calculate the distance between the spacecraft and the body
        const vec3 delta_pos = vec3_sub(gb->bodies[i].pos, craft->pos);
        const double r_squared = vec3_mag_sq(delta_pos);
        if (r_squared < closest_r_squared) {
            closest_r_squared = r_squared;
            closest_planet_id = i;
        }
    }
    craft->closest_r_squared = closest_r_squared;
    craft->closest_planet_id = closest_planet_id;
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

void craft_consumeFuel(spacecraft_t* craft, const double dt) {
    if (craft->engine_on && craft->fuel_mass > 0) {
        double fuel_consumed = craft->mass_flow_rate * craft->throttle * dt;

        if (fuel_consumed > craft->fuel_mass) {
            fuel_consumed = craft->fuel_mass;
            craft->engine_on = false;//
        }

        craft->fuel_mass -= fuel_consumed;
        craft->current_total_mass = craft->dry_mass + craft->fuel_mass;
    }
}

// adds a spacecraft to the spacecraft array
void craft_addSpacecraft(spacecraft_properties_t* gs, const char* name,
                        const vec3 pos, const vec3 vel,
                        const double dry_mass, const double fuel_mass, const double thrust,
                        const double specific_impulse, const double mass_flow_rate,
                        const double attitude_angle, const double moment_of_inertia,
                        const double nozzle_gimbal_range,
                        const burn_properties_t* burns, const int num_burns) {

    // grow capacity if needed
    if (gs->count >= gs->capacity) {
        const int new_capacity = gs->capacity == 0 ? 4 : gs->capacity * 2;
        spacecraft_t* temp = (spacecraft_t*)realloc(gs->spacecraft, new_capacity * sizeof(spacecraft_t));
        if (temp == NULL) {
            displayError("ERROR", "Failed to allocate memory for spacecraft");
            return;
        }
        gs->spacecraft = temp;
        gs->capacity = new_capacity;
    }

    const int idx = gs->count;
    spacecraft_t* craft = &gs->spacecraft[idx];

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
    craft->SOI_planet_id = 0;
    craft->closest_r_squared = INFINITY;
    craft->closest_planet_id = 0;

    craft->apoapsis = 0.0;
    craft->periapsis = 0.0;
    craft->semi_major_axis = 0.0;
    craft->eccentricity = 0.0;

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

    gs->count++;
}
