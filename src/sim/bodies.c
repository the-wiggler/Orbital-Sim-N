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
void body_calculateGravForce(sim_properties_t* sim, const int i, const int j) {
    body_t* bi = &sim->gb.bodies[i];
    body_t* bj = &sim->gb.bodies[j];

    // calculate the distance between the two bodies
    const vec3 delta_pos = vec3_sub(bj->pos, bi->pos);
    const double r_squared = vec3_mag_sq(delta_pos);

    // planet collision logic -- checks if planets are too close
    const double radius_squared = bi->radius * bi->radius;
    if (r_squared < radius_squared) {
        sim->wp.sim_running = false;
        sim->wp.reset_sim = true;
        char err_txt[128];
        snprintf(err_txt, sizeof(err_txt), "Warning: %s has collided with %s\n\nResetting Simulation...", bi->name, bj->name);
        displayError("PLANET COLLISION", err_txt);
        return;
    }

    // force = (G * m1 * m2) * delta / r^3
    const double r = sqrt(r_squared);
    const double r_cubed = r_squared * r;
    const double force_factor = (G * bi->mass * bj->mass) / r_cubed;

    const vec3 force = vec3_scale(delta_pos, force_factor);

    // applies force to both bodies (Newton's third law)
    bi->force = vec3_add(bi->force, force);
    bj->force = vec3_sub(bj->force, force);
}

// calculates changes of velocity and position based on force values
// uses velocity verlet integration
void body_updateMotion(body_t* body, const double dt) {
    // calculate the current acceleration from the force on the object
    body->acc = vec3_scale(body->force, 1.0 / body->mass);

    // update position using current velocity and acceleration
    const vec3 vel_term = vec3_scale(body->vel, dt);
    const vec3 acc_term = vec3_scale(body->acc, 0.5 * dt * dt);
    body->pos = vec3_add(body->pos, vec3_add(vel_term, acc_term));

    // update velocity using average of current and previous acceleration
    const vec3 avg_acc = vec3_scale(vec3_add(body->acc, body->acc_prev), 0.5);
    body->vel = vec3_add(body->vel, vec3_scale(avg_acc, dt));
    body->vel_mag = vec3_mag(body->vel);

    // store current acceleration for next iteration
    body->acc_prev = body->acc;
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
void body_calculateSOI(const body_properties_t* gb) {
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

// function to add a new body to the system
void body_addOrbitalBody(body_properties_t* gb, const char* name, const double mass,
                         const double radius, const vec3 pos, const vec3 vel) {
    // grow capacity if needed (amortized growth)
    if (gb->count >= gb->capacity) {
        const int new_capacity = gb->capacity == 0 ? 4 : gb->capacity * 2;
        body_t* temp = (body_t*)realloc(gb->bodies, new_capacity * sizeof(body_t));
        if (temp == NULL) {
            displayError("ERROR", "Failed to allocate memory for body");
            return;
        }
        gb->bodies = temp;
        gb->capacity = new_capacity;
    }

    const int idx = gb->count;
    body_t* body = &gb->bodies[idx];

    // allocate and copy name
    body->name = (char*)malloc(strlen(name) + 1);
    if (body->name == NULL) {
        displayError("ERROR", "Error: Failed to allocate memory for body name\n");
        return;
    }
#ifdef WIN32
    strcpy_s(body->name, strlen(name) + 1, name);
#else
    strcpy(body->name, name);
#endif

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
