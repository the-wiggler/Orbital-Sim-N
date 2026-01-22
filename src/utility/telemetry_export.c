//
// Created by Java on 12/23/25.
//

#include "telemetry_export.h"
#include <stdio.h>

void writeCSVHeader(FILE* file) {
    fprintf(file, "timestamp,body_index,pos_x,pos_y,vel_x,vel_y,acc_x,acc_y,force_x,force_y\n");
}

void exportTelemetryCSV(const binary_filenames_t filenames, const sim_properties_t sim) {
    const body_properties_t* gb = &sim.gb;
    const window_params_t* wp = &sim.wp;

    // write body position data to the CSV file
    for (int i = 0; i < gb->count; i++) {
        const body_t* body = &gb->bodies[i];

        fprintf(filenames.global_data_FILE, "%f,%d,%f,%f,%f,%f,%f,%f,%f,%f\n",
                wp->sim_time,
                i,
                body->pos.x,
                body->pos.y,
                body->vel.x,
                body->vel.y,
                body->acc.x,
                body->acc.y,
                body->force.x,
                body->force.y);
    }
}
