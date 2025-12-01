#include "config.h"
#include "sim_calculations.h"
#include "sdl_elements.h"
#include <stdio.h>
#include <stdlib.h>
#include <cjson/cJSON.h>

// json handling logic for reading json files
void readBodyJSON(const char* FILENAME, body_properties_t* gb) {
    FILE *fp = fopen(FILENAME, "r");
    if (fp == NULL) {
        displayError("ERROR", "Error: Could not open Body Properties JSON File");
        return;
    }

    // read entire file into buffer
    fseek(fp, 0, SEEK_END);
    const long file_size = ftell(fp);
    fseek(fp, 0, SEEK_SET);

    char* json_buffer = (char*)malloc(file_size + 1);
    if (json_buffer == NULL) {
        displayError("ERROR", "Error: Could not allocate memory for JSON buffer\n");
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
        displayError("ERROR", "Error: Failed to parse JSON from Body JSON");
        return;
    }

    // get bodies array
    const cJSON* bodies = cJSON_GetObjectItemCaseSensitive(json, "bodies");
    if (bodies == NULL || !cJSON_IsArray(bodies)) {
        displayError("ERROR", "Error: 'bodies' array not found or invalid in bodies JSON");
        cJSON_Delete(json);
        return;
    }

    // iterate through bodies
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

    cJSON_Delete(json);
}

void readSpacecraftJSON(const char* FILENAME, spacecraft_properties_t* sc) {
    FILE *fp = fopen(FILENAME, "r");
    if (fp == NULL) {
        displayError("ERROR", "Error: Could not open spacecraft JSON");
        return;
    }

    // read entire file into buffer
    fseek(fp, 0, SEEK_END);
    const long file_size = ftell(fp);
    fseek(fp, 0, SEEK_SET);

    char* json_buffer = (char*)malloc(file_size + 1);
    if (json_buffer == NULL) {
        displayError("ERROR", "Error: Could not allocate memory for body JSON buffer");
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
        displayError("ERROR", "Error: Failed to parse spacecraft file JSON");
        return;
    }

    // get spacecraft array
    const cJSON* spacecraft = cJSON_GetObjectItemCaseSensitive(json, "spacecraft");
    if (spacecraft == NULL || !cJSON_IsArray(spacecraft)) {
        displayError("ERROR", "Error: 'spacecraft' array not found or invalid");
        cJSON_Delete(json);
        return;
    }

    // iterate through spacecraft
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
        cJSON* burn_start_time_item = cJSON_GetObjectItemCaseSensitive(craft, "burn_start_time");
        cJSON* burn_duration_item = cJSON_GetObjectItemCaseSensitive(craft, "burn_duration");
        cJSON* burn_heading_item = cJSON_GetObjectItemCaseSensitive(craft, "burn_heading");
        cJSON* burn_throttle_item = cJSON_GetObjectItemCaseSensitive(craft, "burn_throttle");

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
                            burn_start_time_item->valuedouble,
                            burn_duration_item->valuedouble,
                            burn_heading_item->valuedouble,
                            burn_throttle_item->valuedouble);
    }

    cJSON_Delete(json);
}