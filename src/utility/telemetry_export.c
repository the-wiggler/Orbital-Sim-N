//
// Created by Java on 12/23/25.
//

#include "telemetry_export.h"

void exportTelemetryBinary(binary_filenames_t filenames, window_params_t wp, body_properties_t gb, spacecraft_properties_t sc) {
    // write body position data to the .bin file if enabled
    for (int i = 0; i < gb.count; i++) {
        body_pos_data bpd;
        bpd.timestamp = wp.sim_time;
        bpd.body_index = i;
        bpd.pos_data_x = gb.pos_x[i];
        bpd.pos_data_y = gb.pos_y[i];

        fwrite(&bpd, sizeof(bpd), 1, filenames.body_pos_FILE);
    }
}