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
        global_data_t gd;
        gd.timestamp = wp->sim_time;
        gd.body_index = i;
        gd.pos_data_x = gb->pos_x[i];
        gd.pos_data_y = gb->pos_y[i];
        gd.vel_data_x = gb->vel_x[i];
        gd.vel_data_y = gb->vel_y[i];
        gd.acc_data_x = gb->acc_x[i];
        gd.acc_data_y = gb->acc_y[i];
        gd.force_data_x = gb->force_x[i];
        gd.force_data_y = gb->force_y[i];

        fwrite(&gd, sizeof(gd), 1, filenames.global_data_FILE);
    }
}