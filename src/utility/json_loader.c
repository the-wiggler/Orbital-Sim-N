#include "../utility/json_loader.h"

#include <stdio.h>
#include <string.h>
#include <cjson/cJSON.h>

#include "../globals.h"
#include "../sim/bodies.h"
#include "../sim/spacecraft.h"
#include "../types.h"
#include "../math/matrix.h"

void displayError(const char* title, const char* message);

int findBurnTargetID(const body_properties_t* global_bodies, const char* target_name) {
    for (int i = 0; i < global_bodies->count; i++) {
        if (strcmp(global_bodies->bodies[i].name, target_name) == 0) {
            return i;
        }
    }
    displayError("ERROR", "Relative burn target body not found");
    return -1;
}

relative_burn_target_t findRelativeBurnType(const char* input_burn_type) {
    if (strcmp(input_burn_type, "tangent") == 0) {
        return (relative_burn_target_t){.tangent = true};
    }
    if (strcmp(input_burn_type, "normal") == 0) {
        return (relative_burn_target_t){.normal = true};
    }
    if (strcmp(input_burn_type, "absolute") == 0) {
        return (relative_burn_target_t){.absolute = true};
    }

    displayError("ERROR", "Invalid relative burn type");
    return (relative_burn_target_t){0};
}

// finds the body position by name for relative spacecraft positioning
vec3 findBodyPosition(const body_properties_t* global_bodies, const char* body_name) {
    for (int i = 0; i < global_bodies->count; i++) {
        if (strcmp(global_bodies->bodies[i].name, body_name) == 0) {
            return global_bodies->bodies[i].pos;
        }
    }
    displayError("ERROR", "Body not found for relative positioning");
    return vec3_zero();
}

// finds the body velocity by name for relative spacecraft positioning
vec3 findBodyVelocity(const body_properties_t* global_bodies, const char* body_name) {
    for (int i = 0; i < global_bodies->count; i++) {
        if (strcmp(global_bodies->bodies[i].name, body_name) == 0) {
            return global_bodies->bodies[i].vel;
        }
    }
    return vec3_zero();
}

// json handling logic for reading json files
// NOLINTNEXTLINE(readability-function-cognitive-complexity)
void readSimulationJSON(const char* FILENAME, body_properties_t* global_bodies, spacecraft_properties_t* global_spacecraft) {
    #ifdef _WIN32
    FILE *file_ptr;
    fopen_s(&file_ptr, FILENAME, "r");
    #else
    FILE *file_ptr = fopen(FILENAME, "r");
    #endif
    if (file_ptr == NULL) {
        displayError("ERROR", "Error: Could not open simulation JSON file");
        return;
    }

    // read entire file into buffer
    fseek(file_ptr, 0, SEEK_END);
    const long file_size = ftell(file_ptr);
    fseek(file_ptr, 0, SEEK_SET);

    // validate file size
    if (file_size >= JSON_BUFFER_SIZE) {
        char err_txt[128];
        snprintf(err_txt, sizeof(err_txt), "JSON file too large: %ld bytes (max %d)", file_size, JSON_BUFFER_SIZE);
        displayError("ERROR", err_txt);
        fclose(file_ptr);
        return;
    }

    char json_buffer[JSON_BUFFER_SIZE];
    fread(json_buffer, 1, (size_t)file_size, file_ptr);
    // NOLINTNEXTLINE(clang-analyzer-security.ArrayBound) - safe because we check file_size < JSON_BUFFER_SIZE above
    json_buffer[file_size] = '\0';
    fclose(file_ptr);

    // parse json
    cJSON* json = cJSON_Parse(json_buffer);

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
            cJSON* pos_z_item = cJSON_GetObjectItemCaseSensitive(body, "pos_z");
            cJSON* vel_x_item = cJSON_GetObjectItemCaseSensitive(body, "vel_x");
            cJSON* vel_y_item = cJSON_GetObjectItemCaseSensitive(body, "vel_y");
            cJSON* vel_z_item = cJSON_GetObjectItemCaseSensitive(body, "vel_z");
            cJSON* rotational_v_item = cJSON_GetObjectItemCaseSensitive(body, "rotational_v");
            cJSON* attitude_axis_x = cJSON_GetObjectItemCaseSensitive(body, "attitude_axis_x");
            cJSON* attitude_axis_y = cJSON_GetObjectItemCaseSensitive(body, "attitude_axis_y");
            cJSON* attitude_axis_z = cJSON_GetObjectItemCaseSensitive(body, "attitude_axis_z");
            cJSON* attitude_angle = cJSON_GetObjectItemCaseSensitive(body, "attitude_angle");

            vec3 pos = {
                pos_x_item->valuedouble,
                pos_y_item->valuedouble,
                pos_z_item->valuedouble
            };
            vec3 vel = {
                vel_x_item->valuedouble,
                vel_y_item->valuedouble,
                vel_z_item->valuedouble
            };

            body_addOrbitalBody(global_bodies,
                                name_item->valuestring,
                                mass_item->valuedouble,
                                radius_item->valuedouble,
                                pos, vel);

            // set rotational velocity if present in JSON
            body_t* added_body = &global_bodies->bodies[global_bodies->count - 1];
            if (rotational_v_item != NULL) {
                added_body->rotational_v = rotational_v_item->valuedouble;
            }

            // set attitude if present in JSON
            if (attitude_axis_x != NULL && attitude_axis_y != NULL &&
                attitude_axis_z != NULL && attitude_angle != NULL) {
                vec3 axis = {
                    attitude_axis_x->valuedouble,
                    attitude_axis_y->valuedouble,
                    attitude_axis_z->valuedouble
                };
                added_body->attitude = quaternionFromAxisAngle(axis, attitude_angle->valuedouble);
            }
        }
    }

    // calculate SOI for all bodies after they're loaded
    body_calculateSOI(global_bodies);

    // get spacecraft array
    const cJSON* spacecraft = cJSON_GetObjectItemCaseSensitive(json, "spacecraft");
    if (spacecraft != NULL && cJSON_IsArray(spacecraft)) {
        cJSON* craft = NULL;
        cJSON_ArrayForEach(craft, spacecraft) {
            cJSON* name_item = cJSON_GetObjectItemCaseSensitive(craft, "name");
            cJSON* pos_x_item = cJSON_GetObjectItemCaseSensitive(craft, "pos_x");
            cJSON* pos_y_item = cJSON_GetObjectItemCaseSensitive(craft, "pos_y");
            cJSON* pos_z_item = cJSON_GetObjectItemCaseSensitive(craft, "pos_z");
            cJSON* vel_x_item = cJSON_GetObjectItemCaseSensitive(craft, "vel_x");
            cJSON* vel_y_item = cJSON_GetObjectItemCaseSensitive(craft, "vel_y");
            cJSON* vel_z_item = cJSON_GetObjectItemCaseSensitive(craft, "vel_z");
            cJSON* position_relative_to_item = cJSON_GetObjectItemCaseSensitive(craft, "position_relative_to");
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
            burn_properties_t burns[MAX_BURNS_PER_SPACECRAFT] = {0};

            if (burns_array != NULL && cJSON_IsArray(burns_array)) {
                num_burns = cJSON_GetArraySize(burns_array);
                if (num_burns > MAX_BURNS_PER_SPACECRAFT) {
                    char err_txt[128];
                    snprintf(err_txt, sizeof(err_txt), "Warning: Spacecraft has %d burns, exceeding maximum of %d. Truncating.",
                             num_burns, MAX_BURNS_PER_SPACECRAFT);
                    displayError("WARNING", err_txt);
                    num_burns = MAX_BURNS_PER_SPACECRAFT;
                }
                if (num_burns > 0) {
                    cJSON* burn = NULL;
                    int burn_idx = 0;
                    cJSON_ArrayForEach(burn, burns_array) {
                        cJSON* burn_target = cJSON_GetObjectItemCaseSensitive(burn, "burn_target");
                        cJSON* burn_type = cJSON_GetObjectItemCaseSensitive(burn, "burn_type");
                        cJSON* start_time = cJSON_GetObjectItemCaseSensitive(burn, "start_time");
                        cJSON* duration = cJSON_GetObjectItemCaseSensitive(burn, "duration");
                        cJSON* heading = cJSON_GetObjectItemCaseSensitive(burn, "heading");
                        cJSON* throttle = cJSON_GetObjectItemCaseSensitive(burn, "throttle");

                        const char* burn_target_name = burn_target->valuestring;
                        const int burn_target_id = findBurnTargetID(global_bodies, burn_target_name);
                        if (burn_target_id == -1) {
                            char err_txt[64];
                            snprintf(err_txt, sizeof(err_txt), "Burn target %s not found or is invalid", burn_target_name);
                            displayError("ERROR", err_txt);
                            break;
                        }

                        const char* burn_type_name = burn_type->valuestring;
                        const relative_burn_target_t relative_burn_target_type = findRelativeBurnType(burn_type_name);

                        burns[burn_idx].burn_target_id = burn_target_id;
                        burns[burn_idx].relative_burn_target = relative_burn_target_type;
                        burns[burn_idx].burn_start_time = start_time->valuedouble;
                        burns[burn_idx].burn_end_time = start_time->valuedouble + duration->valuedouble;
                        burns[burn_idx].burn_heading = heading->valuedouble;
                        burns[burn_idx].throttle = throttle->valuedouble;
                        burn_idx++;
                    }
                }
            }

            // calculate final position and velocity for the craft
            vec3 final_pos = {
                pos_x_item->valuedouble,
                pos_y_item->valuedouble,
                pos_z_item->valuedouble
            };
            vec3 final_vel = {
                vel_x_item->valuedouble,
                vel_y_item->valuedouble,
                vel_z_item->valuedouble
            };

            if (position_relative_to_item != NULL && cJSON_IsString(position_relative_to_item)) {
                const char* relative_to = position_relative_to_item->valuestring;
                if (strcmp(relative_to, "absolute") != 0) {
                    // if the position is relative to a body, add its position and velocity
                    vec3 body_pos = findBodyPosition(global_bodies, relative_to);
                    vec3 body_vel = findBodyVelocity(global_bodies, relative_to);
                    final_pos = vec3_add(final_pos, body_pos);
                    final_vel = vec3_add(final_vel, body_vel);
                }
            }

            craft_addSpacecraft(global_spacecraft,
                                name_item->valuestring,
                                final_pos, final_vel,
                                dry_mass_item->valuedouble,
                                fuel_mass_item->valuedouble,
                                thrust_item->valuedouble,
                                specific_impulse_item->valuedouble,
                                mass_flow_rate_item->valuedouble,
                                attitude_item->valuedouble,
                                moment_of_inertia_item->valuedouble,
                                nozzle_gimbal_range_item->valuedouble,
                                burns, num_burns);
        }
    }
    // set the initial closest planet on initialization
    for (int i = 0; i < global_spacecraft->count; i++) {
        craft_findClosestPlanet(&global_spacecraft->spacecraft[i], global_bodies);
    }
    for (int i = 0; i < global_bodies->count; i++) {
        body_findClosestPlanet(&global_bodies->bodies[i], global_bodies);
    }

    cJSON_Delete(json);
}
