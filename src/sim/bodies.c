#include "bodies.h"
#include "../globals.h"
#include "../math/matrix.h"
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

void displayError(const char* title, const char* message);

// calculates gravitational force between two bodies and applies it to both
// i is the body that has the force applied to it, whilst j is the body applying force to i
void body_calculateGravForce(sim_properties_t* sim, const int force_recipient, const int force_applier) {
    body_t* recipient_body = &sim->gb.bodies[force_recipient]; // the body recieving the applied force
    body_t* application_body = &sim->gb.bodies[force_applier]; // the body applying the force

    // calculate the distance between the two bodies
    const vec3 delta_pos = vec3_sub(application_body->pos, recipient_body->pos);
    const double r_squared = vec3_mag_sq(delta_pos);

    // planet collision logic -- checks if planets are too close (sum of both radii)
    const double combined_radius = recipient_body->radius + application_body->radius;
    const double combined_radius_squared = combined_radius * combined_radius;
    if (r_squared < combined_radius_squared) {
        sim->wp.sim_running = false;
        sim->wp.reset_sim = true;
        char err_txt[128];
        snprintf(err_txt, sizeof(err_txt), "Warning: %s has collided with %s\n\nResetting Simulation...", recipient_body->name, application_body->name);
        displayError("PLANET COLLISION", err_txt);
        return;
    }

    // force = (G * m1 * m2) * delta / r^3
    const double r = sqrt(r_squared);
    const double r_cubed = r_squared * r;
    const double force_factor = (G * recipient_body->mass * application_body->mass) / r_cubed;

    const vec3 force = vec3_scale(delta_pos, force_factor);

    // applies force to both bodies
    recipient_body->force = vec3_add(recipient_body->force, force);
    application_body->force = vec3_sub(application_body->force, force);

    // closest planet and SOI tracking for recipient body
    if (r_squared < recipient_body->oe.closest_r_squared) {
        recipient_body->oe.closest_r_squared = r_squared;
        recipient_body->oe.closest_planet_id = force_applier;
        if (r <= application_body->SOI_radius) {
            recipient_body->oe.SOI_planet_id = force_applier;
        }
    }

    // closest planet and SOI tracking for applier body
    if (r_squared < application_body->oe.closest_r_squared) {
        application_body->oe.closest_r_squared = r_squared;
        application_body->oe.closest_planet_id = force_recipient;
        if (r <= recipient_body->SOI_radius) {
            application_body->oe.SOI_planet_id = force_recipient;
        }
    }
}

// calculates the kinetic energy of a target body
void body_calculateKineticEnergy(body_t* body) {
    // calculate kinetic energy (0.5mv^2)
    body->kinetic_energy = 0.5 * body->mass * body->vel_mag * body->vel_mag;
}

// updates the rotational attitude of a body based on its rotational velocity
void body_updateRotation(body_t* body, const double dt) {
    if (body->rotational_v != 0.0) {
        // extract the rotation axis from the current attitude
        const vec3 local_z = {0.0, 0.0, 1.0};
        const vec3 world_spin_axis = quaternionRotate(body->attitude, local_z);

        // create a rotation around this world-space axis
        const double rotation_angle = body->rotational_v * dt;
        const quaternion_t delta_rotation = quaternionFromAxisAngle(world_spin_axis, rotation_angle);

        // apply rotation in world frame
        body->attitude = quaternionMul(delta_rotation, body->attitude);
    }
}

// calculates the SOI radius for all bodies
// assumes the first body (index 0) is the central body
// SOI = a * (m/M)^(2/5) where a is semi-major axis, m is body mass, M is parent mass
void body_calculateSOI(body_properties_t* gb) {
    if (gb->count < 2) return;  // need at least 2 bodies

    body_t* central = &gb->bodies[0];
    const double M = central->mass;

    // first body has no SOI
    central->SOI_radius = 0.0;

    for (int i = 1; i < gb->count; i++) {
        body_t* body = &gb->bodies[i];

        // calculate distance from central body
        const vec3 delta = vec3_sub(body->pos, central->pos);
        const double a = vec3_mag(delta);

        // SOI = a * (m/M)^(2/5)
        const double mass_ratio = body->mass / M;
        body->SOI_radius = a * pow(mass_ratio, 0.4);
    }
}

void body_findClosestPlanet(body_t* body, const body_properties_t* gb) {
    double closest_r_squared = INFINITY;
    int closest_planet_id = 0;
    for (int i = 0; i < gb->count; i++) {
        // calculate the distance between the spacecraft and the body
        const vec3 delta_pos = vec3_sub(gb->bodies[i].pos, body->pos);
        const double r_squared = vec3_mag_sq(delta_pos);
        if (r_squared < closest_r_squared) {
            closest_r_squared = r_squared;
            closest_planet_id = i;
        }
    }
    body->oe.closest_r_squared = closest_r_squared;
    body->oe.closest_planet_id = closest_planet_id;
}

// function to add a new body to the system
void body_addOrbitalBody(body_properties_t* gb, const char* name, const double mass,
                         const double radius, const vec3 pos, const vec3 vel) {
    // check bounds
    if (gb->count >= MAX_PLANETS) {
        char err_txt[128];
        snprintf(err_txt, sizeof(err_txt), "Cannot add body '%s': Maximum of %d bodies reached", name, MAX_PLANETS);
        displayError("ERROR", err_txt);
        return;
    }

    const int idx = gb->count;
    body_t* body = &gb->bodies[idx];

    // validate and copy name
    const size_t name_len = strlen(name);
    if (name_len >= MAX_NAME_LENGTH) {
        char err_txt[128];
        snprintf(err_txt, sizeof(err_txt), "Warning: Body name '%s' truncated to %d characters", name, MAX_NAME_LENGTH - 1);
        displayError("WARNING", err_txt);
    }
#ifdef _WIN32
    strncpy_s(body->name, MAX_NAME_LENGTH, name, MAX_NAME_LENGTH - 1);
#else
    strncpy(body->name, name, MAX_NAME_LENGTH - 1);
#endif
    body->name[MAX_NAME_LENGTH - 1] = '\0';

    // initialize properties
    body->mass = mass;
    body->radius = radius;
    body->SOI_radius = 0.0;
    body->pixel_radius = 0.0f;

    body->pos = pos;
    body->vel = vel;
    body->vel_mag = vec3_mag(vel);
    body->acc = vec3_zero();
    body->acc_prev = vec3_zero();
    body->force = vec3_zero();

    body->kinetic_energy = 0.5 * mass * body->vel_mag * body->vel_mag;
    body->rotational_v = 0.0;
    body->attitude = (quaternion_t){1.0, 0.0, 0.0, 0.0};

    gb->count++;
}
