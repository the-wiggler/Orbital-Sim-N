////////////////////////////////////////////////////////////////////////////////////////////////////
// SIMULATION CALCULATION FUNCTIONS
////////////////////////////////////////////////////////////////////////////////////////////////////
#include "calculation_functions.h"

// at the end of each sim loop, this function should be run to calculate the changes in
// the force values based on other parameters. for example, using F to find a based on m.
// b is the body that has the force applied to it, whilst b2 is the body applying force to b
void calculateForce(body_properties_t *b, body_properties_t b2) {
    // calculate the distance between the two bodies
    double delta_pos_x = b2.pos_x - b->pos_x;
    double delta_pos_y = b2.pos_y - b->pos_y;
    b->r_from_body = sqrt(delta_pos_x * delta_pos_x + delta_pos_y * delta_pos_y);
    double r = b->r_from_body;

    // calculate the force that b2 applies on b due to gravitation (F = (GMm) / r)
    double total_force = (G * b->mass * b2.mass) / (r * r);
    b->force_x += total_force * (delta_pos_x / r);
    b->force_y += total_force * (delta_pos_y / r);
}

// this calculates the changes of velocity and position based on the force values from before
void updateMotion(body_properties_t *b, double dt) {
    // calculate the acceleration from the force on the object
    b->acc_x = b->force_x / b->mass;
    b->acc_y = b->force_y / b->mass;

    // update the velocity
    b->vel_x = b->vel_x + b->acc_x * dt;
    b->vel_y = b->vel_y + b->acc_y * dt;
    b->vel   = sqrt(b->vel_x * b->vel_x + b->vel_y * b->vel_y);

    // update the position using the new velocity
    b->pos_x = b->pos_x + b->vel_x * dt;
    b->pos_y = b->pos_y + b->vel_y * dt;
}

// transforms spacial coordinates (for example, in meters) to pixel coordinates
void transformCoordinates(body_properties_t *b, window_params_t wp) {
    b->pixel_coordinates_x = wp.screen_origin_x + (int)(b->pos_x / wp.meters_per_pixel);
    b->pixel_coordinates_y = wp.screen_origin_y - (int)(b->pos_y / wp.meters_per_pixel); // this is negative because the SDL origin is in the top left, so positive y is 'down'
}

// calculates the size (in pixels) that the planet should appear on the screen based on its mass
int calculateVisualRadius(body_properties_t body, window_params_t wp) {
    int r = (int)(body.radius / wp.meters_per_pixel);
    return r;
}

// function to add a new body to the system
void addOrbitalBody(body_properties_t** gb, int* num_bodies, double mass, double x_pos, double y_pos, double x_vel, double y_vel) {

    // reallocate memory for the new body
    *gb = (body_properties_t *)realloc(*gb, (*num_bodies + 1) * sizeof(body_properties_t));

    // initialize the new body at index num_bodies
    (*gb)[*num_bodies].mass = mass;
    (*gb)[*num_bodies].pos_x = x_pos;
    (*gb)[*num_bodies].pos_y = y_pos;
    (*gb)[*num_bodies].vel_x = x_vel;
    (*gb)[*num_bodies].vel_y = y_vel;
    
    // calculate the radius based on mass
    (*gb)[*num_bodies].radius = pow(mass, 0.279f);

    // increment the body count
    (*num_bodies)++;
}

// reset the simulation by removing all bodies from the system
void resetSim(double* sim_time, body_properties_t** gb, int* num_bodies) {
    // reset simulation time to 0
    *sim_time = 0;

    // free all bodies from memory
    if (*gb != NULL) {
        free(*gb);
        *gb = NULL;
    }

    // reset body count to 0
    *num_bodies = 0;
}

// calculate the optimum velocity for an object to orbit a given body based on the orbit radius


////////////////////////////////////////////////////////////////////////////////////////////////////
// MAIN CALCULATION LOOP
////////////////////////////////////////////////////////////////////////////////////////////////////
void runCalculations(body_properties_t** gb, window_params_t* wp, int num_bodies) {
    if (wp->sim_running) {
        // calculate forces between all body pairs
        if (gb != NULL && *gb != NULL) {
            for (int i = 0; i < num_bodies; i++) {
                (*gb)[i].force_x = 0;
                (*gb)[i].force_y = 0;
                for (int j = 0; j < num_bodies; j++) {
                    if (i != j) {
                        calculateForce(&(*gb)[i], (*gb)[j]);
                    }
                }
            }

            // update the motion for each body and draw
            for (int i = 0; i < num_bodies; i++) {
                // updates the kinematic properties of each body (velocity, accelertion, position, etc)
                updateMotion(&(*gb)[i], wp->time_step);
                // transform real-space coordinate to pixel coordinates on screen (scaling)
                transformCoordinates(&(*gb)[i], *wp);
            }
        }
        wp->sim_time += wp->time_step;
    }
}