//
// Created by Java on 12/23/25.
//

#ifndef ORBITSIMULATION_TELEMETRY_EXPORT_H
#define ORBITSIMULATION_TELEMETRY_EXPORT_H

#include "../types.h"

void exportTelemetryCSV(binary_filenames_t filenames, sim_properties_t sim);
void writeCSVHeader(FILE* file);

#endif //ORBITSIMULATION_TELEMETRY_EXPORT_H