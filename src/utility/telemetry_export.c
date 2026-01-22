//
// Created by Java on 12/23/25.
//

#include "telemetry_export.h"
#include <stdio.h>

void writeCSVHeader(FILE* file) {
    fprintf(file, "timestamp,body_name,planet_x,planet_y,planet_z,planet_vx,planet_vy,planet_vz,craft_name,craft_x,craft_y,craft_z,craft_vx,craft_vy,craft_vz\n");
}

void exportTelemetryCSV(const binary_filenames_t filenames, const sim_properties_t sim) {
    const body_properties_t* gb = &sim.gb;
    const spacecraft_properties_t* gs = &sim.gs;
    const window_params_t* wp = &sim.wp;

    // we choose to iterate over the max of body count and spacecraft count
    int max_count = gb->count > gs->count ? gb->count : gs->count;

    for (int i = 0; i < max_count; i++) {
        // write timestamp
        fprintf(filenames.global_data_FILE, "%f,", wp->sim_time);

        // write body data
        if (i < gb->count) {
            const body_t* body = &gb->bodies[i];
            fprintf(filenames.global_data_FILE, "%s,%f,%f,%f,%f,%f,%f,",
                    body->name,
                    body->pos.x, body->pos.y, body->pos.z,
                    body->vel.x, body->vel.y, body->vel.z);
        } else {
            fprintf(filenames.global_data_FILE, ",,,,,,,");
        }

        // write craft data
        if (i < gs->count) {
            const spacecraft_t* craft = &gs->spacecraft[i];
            fprintf(filenames.global_data_FILE, "%s,%f,%f,%f,%f,%f,%f\n",
                    craft->name,
                    craft->pos.x, craft->pos.y, craft->pos.z,
                    craft->vel.x, craft->vel.y, craft->vel.z);
        } else {
            fprintf(filenames.global_data_FILE, ",,,,,,\n");
        }
    }
}
