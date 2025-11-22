////////////////////////////////////////////////////////////////////////////////////////////////////
// SIMULATION CALCULATION FUNCTIONS
////////////////////////////////////////////////////////////////////////////////////////////////////
#include "sim_calculations.h"
#include "config.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


////////////////////////////////////////////////////////////////////////////////////////////////////
// ORBITAL BODY CALCULATIONS
////////////////////////////////////////////////////////////////////////////////////////////////////
// at the end of each sim loop, this function should be run to calculate the changes in
// the force values based on other parameters. for example, using F to find a based on m.
// b is the body that has the force applied to it, whilst b2 is the body applying force to b
void body_calculateGravForce(body_properties_t *b, body_properties_t b2) {
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
void body_updateMotion(body_properties_t *b, double dt) {
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
void body_transformCoordinates(body_properties_t *b, window_params_t wp) {
    b->pixel_coordinates_x = wp.screen_origin_x + (int)(b->pos_x / wp.meters_per_pixel);
    b->pixel_coordinates_y = wp.screen_origin_y - (int)(b->pos_y / wp.meters_per_pixel); // this is negative because the SDL origin is in the top left, so positive y is 'down'
}

// calculates the kinetic energy of a specific body
void body_calculateKineticEnergy(body_properties_t *b) {
    // calculate kinetic energy (0.5mv^2)
    b->kinetic_energy = 0.5f * b->mass * b->vel * b->vel;
}

// determine the potential energy of a body relative to another target body
double body_calculatePotentialEnergy(body_properties_t b, const char* target_name, body_properties_t* gb, int num_bodies) {
    // search for the target body by name
    body_properties_t* target_body = NULL;
    for (int i = 0; i < num_bodies; i++) {
        if (strcmp(gb[i].name, target_name) == 0) {
            target_body = &gb[i];
            break;
        }
    }

    // if target body not found, return 0
    if (target_body == NULL) {
        return 0.0;
    }

    // calculate distance between the bodies
    double delta_pos_x = target_body->pos_x - b.pos_x;
    double delta_pos_y = target_body->pos_y - b.pos_y;
    double r = sqrt(delta_pos_x * delta_pos_x + delta_pos_y * delta_pos_y);

    // avoid division by zero
    if (r == 0.0) {
        return 0.0;
    }

    // calculate gravitational potential energy (U = -GMm/r)
    double potential_energy = -(G * b.mass * target_body->mass) / r;

    return potential_energy;
}

// calculates the size (in pixels) that the planet should appear on the screen based on its mass
int body_calculateVisualRadius(body_properties_t* body, window_params_t wp) {
    int r = (int)(body->radius / wp.meters_per_pixel);
    body->pixel_radius = r;
    return r;
}

// function to add a new body to the system
void body_addOrbitalBody(body_properties_t** gb, int* num_bodies, char* name, double mass, double x_pos, double y_pos, double x_vel, double y_vel) {

    // reallocate memory for the new body
    *gb = (body_properties_t *)realloc(*gb, (*num_bodies + 1) * sizeof(body_properties_t));

    // allocate memory for the name and copy it
    (*gb)[*num_bodies].name = (char*)malloc(strlen(name) + 1);
    strcpy((*gb)[*num_bodies].name, name);

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

////////////////////////////////////////////////////////////////////////////////////////////////////
// SPACECRAFT CALCULATIONS
////////////////////////////////////////////////////////////////////////////////////////////////////
// calculates the force applied on a spacecraft by a specific body
// (this should be used before any other force functions)
void craft_calculateGravForce(spacecraft_properties_t* s, body_properties_t b) {
    // calculate the distance between the spacecraft and the body
    double delta_pos_x = b.pos_x - s->pos_x;
    double delta_pos_y = b.pos_y - s->pos_y;
    double r = sqrt(delta_pos_x * delta_pos_x + delta_pos_y * delta_pos_y);

    // calculate the ship mass with the current amount of fuel in it
    s->current_total_mass = s->fuel_mass + s->dry_mass;

    // calculate the force that the body applies on the craft due to gravitation (F = (GMm) / r)
    double total_force = (G * s->current_total_mass * b.mass) / (r * r);
    s->force_x += total_force * (delta_pos_x / r);
    s->force_y += total_force * (delta_pos_y / r);
}

// updates the force 
void craft_applyThrust(spacecraft_properties_t* s) {
    if (s->engine_on && s->fuel_mass > 0) {
        double current_thrust = s->thrust * s->throttle;
        s->force_x += current_thrust * cos(s->heading);
        s->force_y += current_thrust * sin(s->heading);
    }
}

// check and activate burns
void craft_checkBurnSchedule(spacecraft_properties_t* s, double sim_time) {
    double burn_end_time = s->burn_start_time + s->burn_duration;
    
    // check if within the burn window
    if (sim_time >= s->burn_start_time && sim_time < burn_end_time && s->fuel_mass > 0) {
        s->engine_on = true;
        s->throttle = s->burn_throttle;
        s->heading = s->burn_heading;
    } else if (sim_time >= burn_end_time) {
        s->engine_on = false;
        s->throttle = 0.0;
    }
}


void craft_consumeFuel(spacecraft_properties_t* s, double dt) {
    if (s->engine_on && s->fuel_mass > 0) {
        double fuel_consumed = s->mass_flow_rate * s->throttle * dt;
        
        if (fuel_consumed > s->fuel_mass) {
            fuel_consumed = s->fuel_mass;
            s->engine_on = false;
        }
        
        s->fuel_mass -= fuel_consumed;
        s->current_total_mass = s->dry_mass + s->fuel_mass;
    }
}


// updates the motion of the spacecraft based on the force currently applied to it
void craft_updateMotion(spacecraft_properties_t* s, double dt) {

    // calculate the acceleration from the force on the object
    s->acc_x = s->force_x / s->current_total_mass;
    s->acc_y = s->force_y / s->current_total_mass;

    // update the velocity
    s->vel_x = s->vel_x + s->acc_x * dt;
    s->vel_y = s->vel_y + s->acc_y * dt;
    s->vel   = sqrt(s->vel_x * s->vel_x + s->vel_y * s->vel_y);

    // update the position using the new velocity
    s->pos_x = s->pos_x + s->vel_x * dt;
    s->pos_y = s->pos_y + s->vel_y * dt;
}

// adds a spacecraft to the spacecraft array
void craft_addSpacecraft(spacecraft_properties_t** sc, int* num_craft, char* name,
                        double x_pos, double y_pos, double x_vel, double y_vel,
                        double dry_mass, double fuel_mass, double thrust,
                        double specific_impulse, double mass_flow_rate,
                        double burn_start_time, double burn_duration, 
                        double burn_heading, double burn_throttle) {
    // reallocate memory for the new spacecraft
    *sc = (spacecraft_properties_t *)realloc(*sc, (*num_craft + 1) * sizeof(spacecraft_properties_t));

    // allocate memory for the name and copy it
    (*sc)[*num_craft].name = (char*)malloc(strlen(name) + 1);
    strcpy((*sc)[*num_craft].name, name);

    // initialize the new craft at index num_craft
    (*sc)[*num_craft].pos_x = x_pos;
    (*sc)[*num_craft].pos_y = y_pos;
    (*sc)[*num_craft].vel_x = x_vel;
    (*sc)[*num_craft].vel_y = y_vel;
    (*sc)[*num_craft].vel = sqrt(x_vel * x_vel + y_vel * y_vel);

    // initialize acceleration to zero
    (*sc)[*num_craft].acc_x = 0.0;
    (*sc)[*num_craft].acc_y = 0.0;

    // initialize forces to zero
    (*sc)[*num_craft].force_x = 0.0;
    (*sc)[*num_craft].force_y = 0.0;

    // initialize pixel coordinates
    (*sc)[*num_craft].pixel_coordinates_x = 0;
    (*sc)[*num_craft].pixel_coordinates_y = 0;

    // initialize spacecraft-specific properties
    (*sc)[*num_craft].heading = 0.0;
    (*sc)[*num_craft].dry_mass = dry_mass;
    (*sc)[*num_craft].fuel_mass = fuel_mass;
    (*sc)[*num_craft].current_total_mass = dry_mass + fuel_mass;
    (*sc)[*num_craft].mass_flow_rate = mass_flow_rate;
    (*sc)[*num_craft].thrust = thrust;
    (*sc)[*num_craft].specific_impulse = specific_impulse;
    (*sc)[*num_craft].throttle = 0.0f;
    (*sc)[*num_craft].engine_on = false;

    // initialize burn schedule parameters
    (*sc)[*num_craft].burn_start_time = burn_start_time;
    (*sc)[*num_craft].burn_duration = burn_duration;
    (*sc)[*num_craft].burn_heading = burn_heading;
    (*sc)[*num_craft].burn_throttle = burn_throttle;

    // increment the craft count
    (*num_craft)++;
}

// transform the craft coordinates in meters to pixel coordinates on the screen
void craft_transformCoordinates(spacecraft_properties_t* s, window_params_t wp) {
    s->pixel_coordinates_x = wp.screen_origin_x + (int)(s->pos_x / wp.meters_per_pixel);
    s->pixel_coordinates_y = wp.screen_origin_y - (int)(s->pos_y / wp.meters_per_pixel); // this is negative because the SDL origin is in the top left, so positive y is 'down'
}


////////////////////////////////////////////////////////////////////////////////////////////////////
// LOGIC FUNCTIONS
////////////////////////////////////////////////////////////////////////////////////////////////////
// reset the simulation by removing all bodies from the system
void resetSim(double* sim_time, body_properties_t** gb, int* num_bodies, spacecraft_properties_t** sc, int* num_craft) {
    // reset simulation time to 0
    *sim_time = 0;

    // free all bodies from memory
    if (*gb != NULL) {
        for (int i = 0; i < *num_bodies; i++) {
            free((*gb)[i].name);
        }
        free(*gb);
        *gb = NULL;
    }

    // reset body count to 0
    *num_bodies = 0;

    // free all craft from memory
    if (*sc != NULL) {
        for (int i = 0; i < *num_craft; i++) {
            free((*sc)[i].name);
        }
        free(*sc);
        *sc = NULL;
    }
    *num_craft = 0;
}

// calculate the optimum velocity for an object to orbit a given body based on the orbit radius

// CSV handling logic for implementing orbital bodies via CSV file
void createCSV(char* FILENAME) {
    FILE *fp = fopen(FILENAME, "w");
    fprintf(fp, "Planet Name,mass,pos_x,pos_y,vel_x,vel_y\n");
    fclose(fp);
}

void readCSV(char* FILENAME, body_properties_t** gb, int* num_bodies) {
    FILE *fp = fopen(FILENAME, "r");

    // check if file was opened successfully
    if (fp == NULL) {
        char error_message[512];
        snprintf(error_message, sizeof(error_message), "Failed to open file: %s\n\nMake sure the file exists and is accessible.", FILENAME);
        SDL_ShowSimpleMessageBox(
            SDL_MESSAGEBOX_ERROR,
            "CSV Load Error",
            error_message,
            NULL
        );
        return;
    }

    char buffer[1024];
    int line_count = 0;

    // read the file line by line
    while (fgets(buffer, sizeof(buffer), fp)) {
        line_count++;

        // skip the header line
        if (line_count == 1) {
            continue;
        }

        // parse the CSV line
        char planet_name[256];
        double mass, pos_x, pos_y, vel_x, vel_y;

        // use sscanf to parse the comma-separated values
        int fields_read = sscanf(buffer, "%255[^,],%lf,%lf,%lf,%lf,%lf",
                                 planet_name, &mass, &pos_x, &pos_y, &vel_x, &vel_y);

        // check if all fields were successfully read
        if (fields_read == 6) {
            // add the orbital body to the system
            body_addOrbitalBody(gb, num_bodies, planet_name, mass, pos_x, pos_y, vel_x, vel_y);
        }
    }
    fclose(fp);
}

void readSpacecraftCSV(char* FILENAME, spacecraft_properties_t** sc, int* num_craft) {
    FILE *fp = fopen(FILENAME, "r");

    // check if file was opened successfully
    if (fp == NULL) {
        char error_message[512];
        snprintf(error_message, sizeof(error_message), "Failed to open file: %s\n\nMake sure the file exists and is accessible.", FILENAME);
        SDL_ShowSimpleMessageBox(
            SDL_MESSAGEBOX_ERROR,
            "CSV Load Error",
            error_message,
            NULL
        );
        return;
    }

    char buffer[1024];
    int line_count = 0;

    // read the file line by line
    while (fgets(buffer, sizeof(buffer), fp)) {
        line_count++;

        // skip the header line
        if (line_count == 1) {
            continue;
        }

        // parse the CSV line
        char craft_name[256];
        double pos_x, pos_y, vel_x, vel_y;
        double dry_mass, fuel_mass, thrust, specific_impulse, mass_flow_rate;
        double burn_start_time, burn_duration, burn_heading, burn_throttle;

        int fields_read = sscanf(buffer, "%255[^,],%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf",
                                craft_name, &pos_x, &pos_y, &vel_x, &vel_y,
                                &dry_mass, &fuel_mass, &thrust, &specific_impulse, &mass_flow_rate,
                                &burn_start_time, &burn_duration, &burn_heading, &burn_throttle);

        if (fields_read == 14) {  // Now expecting 14 fields
            craft_addSpacecraft(sc, num_craft, craft_name, pos_x, pos_y, vel_x, vel_y,
                            dry_mass, fuel_mass, thrust, specific_impulse, mass_flow_rate,
                            burn_start_time, burn_duration, burn_heading, burn_throttle);
        }

    }
    fclose(fp);
}


///////////////////////////////////////////////////////////////////////////////////////////////////////
// MAIN CALCULATION LOOP - THIS LOOP IS PLACED IN THE MAIN LOOP, AND DOES ALL NECESSARY CALCULATIONS //
///////////////////////////////////////////////////////////////////////////////////////////////////////
void runCalculations(body_properties_t** gb, spacecraft_properties_t** sc, window_params_t* wp, int num_bodies, int num_craft) {
    if (wp->sim_running) {
        ////////////////////////////////////////////////////////////////
        // calculate forces between all body pairs
        ////////////////////////////////////////////////////////////////
        if (gb != NULL && *gb != NULL) {
            for (int i = 0; i < num_bodies; i++) {
                // initialize forces to zero to re-calculate them
                (*gb)[i].force_x = 0;
                (*gb)[i].force_y = 0;
                // loop through every body and add the resultant force to the subject body force vector
                for (int j = 0; j < num_bodies; j++) {
                    // check if the body is not calculating on itself
                    if (i != j) {
                        body_calculateGravForce(&(*gb)[i], (*gb)[j]);
                    }
                }
            }

            // update the motion for each body and draw
            for (int i = 0; i < num_bodies; i++) {
                // updates the kinematic properties of each body (velocity, accelertion, position, etc)
                body_updateMotion(&(*gb)[i], wp->time_step);
                // transform real-space coordinate to pixel coordinates on screen (scaling)
                body_transformCoordinates(&(*gb)[i], *wp);
            }
        }

        ////////////////////////////////////////////////////////////////
        // calculate forces between spacecraft and bodies
        ////////////////////////////////////////////////////////////////
        if (sc != NULL && *sc != NULL && gb != NULL && *gb != NULL) {
            for (int i = 0; i < num_craft; i++) {
                (*sc)[i].force_x = 0;
                (*sc)[i].force_y = 0;

                // check if burn should be active based on simulation time
                craft_checkBurnSchedule(&(*sc)[i], wp->sim_time);

                // loop through all bodies and calculate gravitational forces on spacecraft
                for (int j = 0; j < num_bodies; j++) {
                    craft_calculateGravForce(&(*sc)[i], (*gb)[j]);
                }

                // apply thrust first, then consume fuel
                craft_applyThrust(&(*sc)[i]);
                craft_consumeFuel(&(*sc)[i], wp->time_step);
            }

            // update the motion for each craft and draw
            for (int i = 0; i < num_craft; i++) {
                // updates the kinematic properties of each body (velocity, accelertion, position, etc)
                craft_updateMotion(&(*sc)[i], wp->time_step);
                // transform real-space coordinate to pixel coordinates on screen (scaling)
                craft_transformCoordinates(&(*sc)[i], *wp);
            }
        }

        // increment the time based on the time step
        wp->sim_time += wp->time_step;
    }
}