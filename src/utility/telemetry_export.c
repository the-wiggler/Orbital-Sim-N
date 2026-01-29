//
// Created by Java on 12/23/25.
//

#include "telemetry_export.h"
#include <stdio.h>
#include "../types.h"

void writeCSVHeader(FILE* file) {
    fprintf(file, "timestamp,body_name,planet_x,planet_y,planet_z,planet_vx,planet_vy,planet_vz,craft_name,craft_x,craft_y,craft_z,craft_vx,craft_vy,craft_vz,,system_energy\n");
}

void exportTelemetryCSV(const binary_filenames_t filenames, const sim_properties_t sim) {
    const body_properties_t* global_bodies = &sim.global_bodies;
    const spacecraft_properties_t* global_spacecraft = &sim.global_spacecraft;
    const window_params_t* window_props = &sim.window_params;

    // we choose to iterate over the max of body count and spacecraft count
    int max_count = global_bodies->count > global_spacecraft->count ? global_bodies->count : global_spacecraft->count;

    for (int i = 0; i < max_count; i++) {
        // write timestamp
        fprintf(filenames.global_data_FILE, "%f,", window_props->sim_time);

        // write body data
        if (i < global_bodies->count) {
            const body_t* body = &global_bodies->bodies[i];
            fprintf(filenames.global_data_FILE, "%s,%f,%f,%f,%f,%f,%f,",
                    body->name,
                    body->pos.x, body->pos.y, body->pos.z,
                    body->vel.x, body->vel.y, body->vel.z);
        } else {
            fprintf(filenames.global_data_FILE, ",,,,,,,");
        }

        // write craft data
        if (i < global_spacecraft->count) {
            const spacecraft_t* craft = &global_spacecraft->spacecraft[i];
            fprintf(filenames.global_data_FILE, "%s,%f,%f,%f,%f,%f,%f",
                    craft->name,
                    craft->pos.x, craft->pos.y, craft->pos.z,
                    craft->vel.x, craft->vel.y, craft->vel.z);
        } else {
            fprintf(filenames.global_data_FILE, ",,,,,,");
        }
        // write energies
        fprintf(filenames.global_data_FILE, ",,%f,\n", sim.global_bodies.total_system_energy);
    }
}
