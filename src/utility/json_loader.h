#ifndef JSON_LOADER_H
#define JSON_LOADER_H

#include "../types.h"

void readSimulationJSON(const char* FILENAME, body_properties_t* global_bodies, spacecraft_properties_t* global_spacecraft);
int findSpacecraftID(const spacecraft_properties_t* global_spacecraft, const char* craft_name);

#endif
