#include "../utility/json_loader.h"
#include "../sim/bodies.h"
#include "../sim/spacecraft.h"
#include <stdio.h>
#include <stdlib.h>
#include <cjson/cJSON.h>

void displayError(const char* title, const char* message);

// json handling logic for reading json files
void readSimulationJSON(const char* FILENAME, body_properties_t* gb, spacecraft_properties_t* sc) {
    #ifdef _WIN32
    FILE *fp;
    fopen_s(&fp, FILENAME, "r");
    #else
    FILE *fp = fopen(FILENAME, "r");
    #endif
    if (fp == NULL) {
        displayError("ERROR", "Error: Could not open simulation JSON file");
        return;
    }

    // read entire file into buffer
    fseek(fp, 0, SEEK_END);
    const long file_size = ftell(fp);
    fseek(fp, 0, SEEK_SET);

    char* json_buffer = (char*)malloc(file_size + 1);
    if (json_buffer == NULL) {
        displayError("ERROR", "Error: Could not allocate memory for JSON buffer");
        fclose(fp);
        return;
    }

    fread(json_buffer, 1, file_size, fp);
    json_buffer[file_size] = '\0';
    fclose(fp);

    // parse json
    cJSON* json = cJSON_Parse(json_buffer);
    free(json_buffer);

    if (json == NULL) {
        displayError("ERROR", "Error: Failed to parse simulation JSON");
        return;
    }

    // get bodies array
    const cJSON* bodies = cJSON_GetObjectItemCaseSensitive(json, "bodies");
    if (bodies != NULL && cJSON_IsArray(bodies)) {
        const cJSON* body = NULL;
        cJSON_ArrayForEach(body, bodies) {
            cJSON* name_item = cJSON_GetObjectItemCaseSensitive(body, "name");
            cJSON* mass_item = cJSON_GetObjectItemCaseSensitive(body, "mass");
            cJSON* radius_item = cJSON_GetObjectItemCaseSensitive(body, "radius");
            cJSON* pos_x_item = cJSON_GetObjectItemCaseSensitive(body, "pos_x");
            cJSON* pos_y_item = cJSON_GetObjectItemCaseSensitive(body, "pos_y");
            cJSON* vel_x_item = cJSON_GetObjectItemCaseSensitive(body, "vel_x");
            cJSON* vel_y_item = cJSON_GetObjectItemCaseSensitive(body, "vel_y");

            body_addOrbitalBody(gb,
                                name_item->valuestring,
                                mass_item->valuedouble,
                                radius_item->valuedouble,
                                pos_x_item->valuedouble,
                                pos_y_item->valuedouble,
                                vel_x_item->valuedouble,
                                vel_y_item->valuedouble);
        }
    }

    // get spacecraft array
    const cJSON* spacecraft = cJSON_GetObjectItemCaseSensitive(json, "spacecraft");
    if (spacecraft != NULL && cJSON_IsArray(spacecraft)) {
        cJSON* craft = NULL;
        cJSON_ArrayForEach(craft, spacecraft) {
            cJSON* name_item = cJSON_GetObjectItemCaseSensitive(craft, "name");
            cJSON* pos_x_item = cJSON_GetObjectItemCaseSensitive(craft, "pos_x");
            cJSON* pos_y_item = cJSON_GetObjectItemCaseSensitive(craft, "pos_y");
            cJSON* vel_x_item = cJSON_GetObjectItemCaseSensitive(craft, "vel_x");
            cJSON* vel_y_item = cJSON_GetObjectItemCaseSensitive(craft, "vel_y");
            cJSON* dry_mass_item = cJSON_GetObjectItemCaseSensitive(craft, "dry_mass");
            cJSON* fuel_mass_item = cJSON_GetObjectItemCaseSensitive(craft, "fuel_mass");
            cJSON* thrust_item = cJSON_GetObjectItemCaseSensitive(craft, "thrust");
            cJSON* specific_impulse_item = cJSON_GetObjectItemCaseSensitive(craft, "specific_impulse");
            cJSON* mass_flow_rate_item = cJSON_GetObjectItemCaseSensitive(craft, "mass_flow_rate");
            cJSON* attitude_item = cJSON_GetObjectItemCaseSensitive(craft, "attitude");
            cJSON* moment_of_inertia_item = cJSON_GetObjectItemCaseSensitive(craft, "moment_of_inertia");
            cJSON* nozzle_gimbal_range_item = cJSON_GetObjectItemCaseSensitive(craft, "nozzle_gimbal_range");

            // parse burns array
            cJSON* burns_array = cJSON_GetObjectItemCaseSensitive(craft, "burns");
            int num_burns = 0;
            burn_properties_t* burns = NULL;

            if (burns_array != NULL && cJSON_IsArray(burns_array)) {
                num_burns = cJSON_GetArraySize(burns_array);
                if (num_burns > 0) {
                    burns = (burn_properties_t*)malloc(num_burns * sizeof(burn_properties_t));
                    cJSON* burn = NULL;
                    int burn_idx = 0;
                    cJSON_ArrayForEach(burn, burns_array) {
                        cJSON* start_time = cJSON_GetObjectItemCaseSensitive(burn, "start_time");
                        cJSON* duration = cJSON_GetObjectItemCaseSensitive(burn, "duration");
                        cJSON* heading = cJSON_GetObjectItemCaseSensitive(burn, "heading");
                        cJSON* throttle = cJSON_GetObjectItemCaseSensitive(burn, "throttle");

                        burns[burn_idx].burn_start_time = start_time->valuedouble;
                        burns[burn_idx].burn_end_time = start_time->valuedouble + duration->valuedouble;
                        burns[burn_idx].burn_heading = heading->valuedouble;
                        burns[burn_idx].throttle = throttle->valuedouble;
                        burn_idx++;
                    }
                }
            }

            craft_addSpacecraft(sc,
                                name_item->valuestring,
                                pos_x_item->valuedouble,
                                pos_y_item->valuedouble,
                                vel_x_item->valuedouble,
                                vel_y_item->valuedouble,
                                dry_mass_item->valuedouble,
                                fuel_mass_item->valuedouble,
                                thrust_item->valuedouble,
                                specific_impulse_item->valuedouble,
                                mass_flow_rate_item->valuedouble,
                                attitude_item->valuedouble,
                                moment_of_inertia_item->valuedouble,
                                nozzle_gimbal_range_item->valuedouble,
                                burns, num_burns);

            // free burns array after adding spacecraft
            if (burns != NULL) {
                free(burns);
            }
        }
    }

    cJSON_Delete(json);
}
