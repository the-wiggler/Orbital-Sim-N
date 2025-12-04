#include "config.h"
#include "sim_calculations.h"
#include "sdl_elements.h"
#include <stdio.h>
#include <stdlib.h>
#include <cjson/cJSON.h>

// json handling logic for reading json files
void readBodyJSON(const char* FILENAME, body_properties_t* gb) {
    #ifdef _WIN32
        FILE *fp;
        fopen_s(&fp, FILENAME, "r");
    #else
        FILE *fp = fopen(FILENAME, "r");
    #endif
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
    #ifdef _WIN32
    FILE *fp;
    fopen_s(&fp, FILENAME, "r");
    #else
    FILE *fp = fopen(FILENAME, "r");
    #endif
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
    }

    cJSON_Delete(json);
}


// cleanup for main
void cleanup(void* args) {
    const cleanup_args* c = (cleanup_args*)args;

    // free all bodies
    for (int i = 0; i < c->gb->count; i++) {
        free(c->gb->names[i]);
    }

    // free cache arrays for each body
    if (c->gb->cached_body_coords_x != NULL) {
        for (int i = 0; i < c->gb->count; i++) {
            free(c->gb->cached_body_coords_x[i]);
        }
        free(c->gb->cached_body_coords_x);
    }
    if (c->gb->cached_body_coords_y != NULL) {
        for (int i = 0; i < c->gb->count; i++) {
            free(c->gb->cached_body_coords_y[i]);
        }
        free(c->gb->cached_body_coords_y);
    }

    free(c->gb->names);
    free(c->gb->mass);
    free(c->gb->radius);
    free(c->gb->pixel_radius);
    free(c->gb->pos_x);
    free(c->gb->pos_y);
    free(c->gb->pixel_coordinates_x);
    free(c->gb->pixel_coordinates_y);
    free(c->gb->vel_x);
    free(c->gb->vel_y);
    free(c->gb->vel);
    free(c->gb->acc_x);
    free(c->gb->acc_y);
    free(c->gb->acc_x_prev);
    free(c->gb->acc_y_prev);
    free(c->gb->force_x);
    free(c->gb->force_y);
    free(c->gb->kinetic_energy);

    // free all spacecraft
    for (int i = 0; i < c->sc->count; i++) {
        free(c->sc->names[i]);
        free(c->sc->burn_properties[i]);
    }
    free(c->sc->names);
    free(c->sc->current_total_mass);
    free(c->sc->dry_mass);
    free(c->sc->fuel_mass);
    free(c->sc->pos_x);
    free(c->sc->pos_y);
    free(c->sc->pixel_coordinates_x);
    free(c->sc->pixel_coordinates_y);
    free(c->sc->attitude);
    free(c->sc->vel_x);
    free(c->sc->vel_y);
    free(c->sc->vel);
    free(c->sc->rotational_v);
    free(c->sc->momentum);
    free(c->sc->acc_x);
    free(c->sc->acc_y);
    free(c->sc->acc_x_prev);
    free(c->sc->acc_y_prev);
    free(c->sc->rotational_a);
    free(c->sc->moment_of_inertia);
    free(c->sc->grav_force_x);
    free(c->sc->grav_force_y);
    free(c->sc->torque);
    free(c->sc->thrust);
    free(c->sc->mass_flow_rate);
    free(c->sc->specific_impulse);
    free(c->sc->throttle);
    free(c->sc->nozzle_gimbal_range);
    free(c->sc->nozzle_velocity);
    free(c->sc->engine_on);
    free(c->sc->num_burns);
    free(c->sc->burn_properties);
}