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

// calculates the kinetic energy of a specific body
void calculateKineticEnergy(body_properties_t *b) {
    // calculate kinetic energy (0.5mv^2)
    b->kinetic_energy = 0.5f * b->mass * b->vel * b->vel;
}

// determine the potential energy of a body relative to another target body
double calculatePotentialEnergy(body_properties_t b, const char* target_name, body_properties_t* gb, int num_bodies) {
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
int calculateVisualRadius(body_properties_t* body, window_params_t wp) {
    int r = (int)(body->radius / wp.meters_per_pixel);
    body->pixel_radius = r;
    return r;
}

// function to add a new body to the system
void addOrbitalBody(body_properties_t** gb, int* num_bodies, char* name, double mass, double x_pos, double y_pos, double x_vel, double y_vel) {

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

// reset the simulation by removing all bodies from the system
void resetSim(double* sim_time, body_properties_t** gb, int* num_bodies) {
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
            addOrbitalBody(gb, num_bodies, planet_name, mass, pos_x, pos_y, vel_x, vel_y);
        }
    }
    fclose(fp);
}

// create a data file
void storeOrbitalData(body_properties_t* gb) {

}

// export orbital energy data between all bodies to a CSV file
// This function appends data to the file, creating it with headers if it doesn't exist
void exportOrbitalEnergyData(char* FILENAME, body_properties_t* gb, int num_bodies, double sim_time) {
    // check if file exists to determine if we need to write headers
    bool file_exists = false;
    FILE *fp_check = fopen(FILENAME, "r");
    if (fp_check != NULL) {
        file_exists = true;
        fclose(fp_check);
    }

    // open file in append mode
    FILE *fp = fopen(FILENAME, "a");

    if (fp == NULL) {
        char error_message[512];
        snprintf(error_message, sizeof(error_message), "Failed to create file: %s\n\nCheck file permissions.", FILENAME);
        SDL_ShowSimpleMessageBox(
            SDL_MESSAGEBOX_ERROR,
            "CSV Export Error",
            error_message,
            NULL
        );
        return;
    }

    // write CSV header only if file is new
    if (!file_exists) {
        fprintf(fp, "Body1,Body2,Body1_Mass,Body2_Mass,Distance,Body1_KineticEnergy,Body2_KineticEnergy,PotentialEnergy,TotalOrbitalEnergy,SimTime\n");
    }

    // calculate kinetic energy for all bodies first
    for (int i = 0; i < num_bodies; i++) {
        calculateKineticEnergy(&gb[i]);
    }

    // iterate through all unique pairs of bodies
    for (int i = 0; i < num_bodies; i++) {
        for (int j = i + 1; j < num_bodies; j++) {
            // calculate distance between bodies
            double delta_pos_x = gb[j].pos_x - gb[i].pos_x;
            double delta_pos_y = gb[j].pos_y - gb[i].pos_y;
            double distance = sqrt(delta_pos_x * delta_pos_x + delta_pos_y * delta_pos_y);

            // calculate potential energy between the pair (U = -GMm/r)
            double potential_energy = 0.0;
            if (distance > 0.0) {
                potential_energy = -(G * gb[i].mass * gb[j].mass) / distance;
            }

            // total orbital energy = KE1 + KE2 + PE
            double total_energy = gb[i].kinetic_energy + gb[j].kinetic_energy + potential_energy;

            // write data row
            fprintf(fp, "%s,%s,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e\n",
                    gb[i].name,
                    gb[j].name,
                    gb[i].mass,
                    gb[j].mass,
                    distance,
                    gb[i].kinetic_energy,
                    gb[j].kinetic_energy,
                    potential_energy,
                    total_energy,
                    sim_time);
        }
    }

    fclose(fp);
}

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

            // log orbital energy data every iteration
            //exportOrbitalEnergyData("orbital_energy.csv", *gb, num_bodies, wp->sim_time);
        }
        wp->sim_time += wp->time_step;
    }
}