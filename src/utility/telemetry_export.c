//
// Created by Java on 12/23/25.
//

#include "telemetry_export.h"
#include <stdio.h>

void exportTelemetryBinary(binary_filenames_t filenames, const sim_properties_t* sim) {
    const body_properties_t* gb = &sim->gb;
    const window_params_t* wp = &sim->wp;

    // write body position data to the .bin file if enabled
    for (int i = 0; i < gb->count; i++) {
        body_t* body = &gb->bodies[i];

        global_data_t gd;
        gd.timestamp = wp->sim_time;
        gd.body_index = i;
        gd.pos_data_x = body->pos.x;
        gd.pos_data_y = body->pos.y;
        gd.vel_data_x = body->vel.x;
        gd.vel_data_y = body->vel.y;
        gd.acc_data_x = body->acc.x;
        gd.acc_data_y = body->acc.y;
        gd.force_data_x = body->force.x;
        gd.force_data_y = body->force.y;

        fwrite(&gd, sizeof(gd), 1, filenames.global_data_FILE);
    }
}
