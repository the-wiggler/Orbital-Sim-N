//
// Created by java on 12/11/2025.
//
#include "../globals.h"
#include "../types.h"
#include "craft_view.h"
#include <math.h>
#include <stdio.h>

#include "renderer.h"

// calculates the distance and heading of a craft to another body
typedef struct {
    double heading;
    double distance_x, distance_y, distance;
    double rel_v_x, rel_v_y, rel_v;
} rel_body_vector_t;
rel_body_vector_t craft_calcVectorToBody(const spacecraft_properties_t* sc, const body_properties_t* gb, const int craft_id, const int body_id) {
    // relative distance
    const double dist_x = gb->pos_x[body_id] - sc->pos_x[craft_id];
    const double dist_y = gb->pos_y[body_id] - sc->pos_y[craft_id];
    const double dist = sqrt(pow(dist_x, 2) + pow(dist_y, 2));

    // relative velocity
    const double rvx = gb->vel_x[body_id] - sc->vel_x[craft_id];
    const double rvy = gb->vel_y[body_id] - sc->vel_y[craft_id];
    const double rv = sqrt(pow(rvx, 2) + pow(rvy, 2));

    // heading (radians)
    const double hd = atan2(-dist_y, dist_x);

    rel_body_vector_t result;
    result.heading = hd;
    result.distance_x = dist_x;
    result.distance_y = dist_y;
    result.distance = dist;
    result.rel_v_x = rvx;
    result.rel_v_y = rvy;
    result.rel_v = rv;
    return result;
}

// draws the shape of the craft in the craft view
void craft_DrawCraft(SDL_Renderer* renderer, const window_params_t* wp, const spacecraft_properties_t* sc, const int craft_id) {
    if (sc != NULL && sc->count > 0 && sc->attitude != NULL) {
        const float center_x = wp->window_size_x / 2.0f;
        const float center_y = wp->window_size_y / 2.0f;

        // arrow shape
        const float arrow_length = 40.0f;
        const float arrow_width = 25.0f;
        const float vertices[3][2] = {
            {0.0f, -arrow_length / 2.0f},
            {-arrow_width / 2.0f, arrow_length / 2.0f},
            {arrow_width / 2.0f, arrow_length / 2.0f}
        };

        // rotate and translate vertices based on heading (radians -- ccw positive)
        SDL_FPoint points[4];
        const float render_attitude =  -(float)sc->attitude[craft_id] + (float)M_PI / 2.0f;
        for (int i = 0; i < 3; i++) {
            const float x = vertices[i][0];
            const float y = vertices[i][1];

            const float rotated_x = x * cosf(render_attitude) - y * sinf(render_attitude);
            const float rotated_y = x * sinf(render_attitude) + y * cosf(render_attitude);

            points[i].x = center_x + rotated_x;
            points[i].y = center_y + rotated_y;
        }
        points[3] = points[0];

        SDL_SetRenderDrawColor(renderer, 255, 255, 255, 255);
        SDL_Vertex triangle_vertices[3];
        for (int i = 0; i < 3; i++) {
            triangle_vertices[i].position = points[i];
            triangle_vertices[i].color.r = 255;
            triangle_vertices[i].color.g = 255;
            triangle_vertices[i].color.b = 255;
            triangle_vertices[i].color.a = 255;
        }

        const int indices[3] = {0, 1, 2};
        SDL_RenderGeometry(renderer, NULL, triangle_vertices, 3, indices, 3);
    }
}

const float vertical_desc_text_buffer = 15;

// draws the arrows that point to planets and their distances relative to the craft
void craft_DrawRelativeArrows(SDL_Renderer* renderer, const window_params_t* wp, const body_properties_t* gb, const spacecraft_properties_t* sc, const int craft_id) {
    if (sc != NULL && sc->count > 0 && gb != NULL && gb->count > 0) {
        const float center_x = wp->window_size_x / 2.0f;
        const float center_y = wp->window_size_y / 2.0f;

        SDL_SetRenderDrawColor(renderer, 255, 255, 255, 255);

        const float arrow_length = 200.0f; // fixed arrow length

        for (int i = 0; i < gb->count; i++) {
            // calculate relative position/heading to the craft
            const rel_body_vector_t body_vector = craft_calcVectorToBody(sc, gb, craft_id, i);

            char dist_text[32];
            snprintf(dist_text, sizeof(dist_text), "%.1f km", (body_vector.distance - gb->radius[i]) / 1000.0f);

            // create a fixed-length arrow pointing in the direction of the body
            const float pointer_ending_x = center_x + arrow_length * cosf((float)body_vector.heading);
            const float pointer_ending_y = center_y + arrow_length * sinf((float)body_vector.heading);

            SDL_RenderLine(renderer, center_x, center_y, pointer_ending_x, pointer_ending_y);
            SDL_WriteText(renderer, g_font_small, gb->names[i], pointer_ending_x, pointer_ending_y, TEXT_COLOR);
            SDL_WriteText(renderer, g_font_small, dist_text, pointer_ending_x, pointer_ending_y + vertical_desc_text_buffer, TEXT_COLOR);
        }
    }
}

// draws the craft's velocity, thrust, etc. arrows
void craft_drawCraftPropertyArrows(SDL_Renderer* renderer, const window_params_t* wp, const body_properties_t* gb, const spacecraft_properties_t* sc, const int craft_id) {
    if (sc != NULL && sc->count > 0 && gb != NULL && gb->count > 0) {
        const float center_x = wp->window_size_x / 2.0f;
        const float center_y = wp->window_size_y / 2.0f;

        const float arrow_length = 100.0f;

        ////////////////////////////////////////////////////
        // absolute velocity arrow relative to space
        ////////////////////////////////////////////////////
        const float relative_abs_v_heading = atan2f(-(float)sc->vel_y[craft_id], (float)sc->vel_x[craft_id]);
        const float abs_v_arrow_end_x = center_x + arrow_length * cosf(relative_abs_v_heading);
        const float abs_v_arrow_end_y = center_y + arrow_length * sinf(relative_abs_v_heading);

        SDL_SetRenderDrawColor(renderer, 255, 0, 0, 255);
        SDL_RenderLine(renderer, center_x, center_y, abs_v_arrow_end_x, abs_v_arrow_end_y);

        char abs_v_text[32];
        snprintf(abs_v_text, sizeof(abs_v_text), "%.2f m/s", sc->vel[craft_id]);
        SDL_WriteText(renderer, g_font_small, "Σv", abs_v_arrow_end_x, abs_v_arrow_end_y, TEXT_COLOR);

        ////////////////////////////////////////////////////
        // net acceleration arrow
        ////////////////////////////////////////////////////
        const float acc_heading = atan2f(-(float)sc->acc_y[craft_id], (float)sc->acc_x[craft_id]);
        const float acc_arrow_x = center_x + arrow_length * cosf(acc_heading);
        const float acc_arrow_y = center_y + arrow_length * sinf(acc_heading);

        SDL_SetRenderDrawColor(renderer, 0, 0, 255, 255);
        SDL_RenderLine(renderer, center_x, center_y, acc_arrow_x, acc_arrow_y);

        char a_text[32];
        snprintf(a_text, sizeof(a_text), "%.2f m/s^2", sc->vel[craft_id]);
        SDL_WriteText(renderer, g_font_small, "Σa", acc_arrow_x, acc_arrow_y, TEXT_COLOR);

    }
}

// craft information
void craft_drawCraftInfo(SDL_Renderer* renderer, const window_params_t* wp, const body_properties_t* gb, const spacecraft_properties_t* sc, const int craft_id) {
    if (sc == NULL || sc->count == 0) return;

    const float padding = wp->window_size_y * 0.03f;
    const float line_height = wp->window_size_y * 0.02f;
    float y_offset = wp->window_size_y * 0.1f;

    // position status bar on the left side
    const float x_pos = padding;

    char text_buffer[128];

    // name
    SDL_WriteText(renderer, g_font_small, sc->names[craft_id], x_pos, y_offset, TEXT_COLOR);
    y_offset += line_height * 1.5f;

    // position
    snprintf(text_buffer, sizeof(text_buffer), "Position:");
    SDL_WriteText(renderer, g_font_small, text_buffer, x_pos, y_offset, TEXT_COLOR);
    y_offset += line_height;

    snprintf(text_buffer, sizeof(text_buffer), "  x: %.2f km", sc->pos_x[craft_id] / 1000.0);
    SDL_WriteText(renderer, g_font_small, text_buffer, x_pos, y_offset, TEXT_COLOR);
    y_offset += line_height;

    snprintf(text_buffer, sizeof(text_buffer), "  y: %.2f km", sc->pos_y[craft_id] / 1000.0);
    SDL_WriteText(renderer, g_font_small, text_buffer, x_pos, y_offset, TEXT_COLOR);
    y_offset += line_height * 1.5f;

    // velocity
    snprintf(text_buffer, sizeof(text_buffer), "Velocity: %.2f m/s", sc->vel[craft_id]);
    SDL_WriteText(renderer, g_font_small, text_buffer, x_pos, y_offset, TEXT_COLOR);
    y_offset += line_height;

    snprintf(text_buffer, sizeof(text_buffer), "  vx: %.2f m/s", sc->vel_x[craft_id]);
    SDL_WriteText(renderer, g_font_small, text_buffer, x_pos, y_offset, TEXT_COLOR);
    y_offset += line_height;

    snprintf(text_buffer, sizeof(text_buffer), "  vy: %.2f m/s", sc->vel_y[craft_id]);
    SDL_WriteText(renderer, g_font_small, text_buffer, x_pos, y_offset, TEXT_COLOR);
    y_offset += line_height * 1.5f;

    // acceleration
    const double acc_mag = sqrt(pow(sc->acc_x[craft_id], 2) + pow(sc->acc_y[craft_id], 2));
    snprintf(text_buffer, sizeof(text_buffer), "Acceleration: %.4f m/s^2", acc_mag);
    SDL_WriteText(renderer, g_font_small, text_buffer, x_pos, y_offset, TEXT_COLOR);
    y_offset += line_height;

    snprintf(text_buffer, sizeof(text_buffer), "  ax: %.4f m/s^2", sc->acc_x[craft_id]);
    SDL_WriteText(renderer, g_font_small, text_buffer, x_pos, y_offset, TEXT_COLOR);
    y_offset += line_height;

    snprintf(text_buffer, sizeof(text_buffer), "  ay: %.4f m/s^2", sc->acc_y[craft_id]);
    SDL_WriteText(renderer, g_font_small, text_buffer, x_pos, y_offset, TEXT_COLOR);
    y_offset += line_height * 1.5f;

    // attitude
    const double attitude_deg = sc->attitude[craft_id] * 180.0 / M_PI;
    snprintf(text_buffer, sizeof(text_buffer), "Attitude: %.1f deg", attitude_deg);
    SDL_WriteText(renderer, g_font_small, text_buffer, x_pos, y_offset, TEXT_COLOR);
    y_offset += line_height;

    snprintf(text_buffer, sizeof(text_buffer), "Rot. vel: %.2f deg/s", sc->rotational_v[craft_id] * 180.0 / M_PI);
    SDL_WriteText(renderer, g_font_small, text_buffer, x_pos, y_offset, TEXT_COLOR);
    y_offset += line_height * 1.5f;

    // mass
    snprintf(text_buffer, sizeof(text_buffer), "Mass: %.0f kg", sc->current_total_mass[craft_id]);
    SDL_WriteText(renderer, g_font_small, text_buffer, x_pos, y_offset, TEXT_COLOR);
    y_offset += line_height;

    snprintf(text_buffer, sizeof(text_buffer), "Fuel: %.0f kg", sc->fuel_mass[craft_id]);
    SDL_WriteText(renderer, g_font_small, text_buffer, x_pos, y_offset, TEXT_COLOR);
    y_offset += line_height * 1.5f;

    // engine status
    snprintf(text_buffer, sizeof(text_buffer), "Engine: %s", sc->engine_on[craft_id] ? "ON" : "OFF");
    SDL_WriteText(renderer, g_font_small, text_buffer, x_pos, y_offset, TEXT_COLOR);
    y_offset += line_height;

    if (sc->engine_on[craft_id]) {
        snprintf(text_buffer, sizeof(text_buffer), "Throttle: %.1f%%", sc->throttle[craft_id] * 100.0);
        SDL_WriteText(renderer, g_font_small, text_buffer, x_pos, y_offset, TEXT_COLOR);
        y_offset += line_height;

        snprintf(text_buffer, sizeof(text_buffer), "Thrust: %.2f N", sc->thrust[craft_id] * sc->throttle[craft_id]);
        SDL_WriteText(renderer, g_font_small, text_buffer, x_pos, y_offset, TEXT_COLOR);
        y_offset += line_height;
    }
    y_offset += line_height * 0.5f;

    // nearest body info
    if (gb != NULL && gb->count > 0) {
        int closest_body_id = 0;
        double min_distance = 0.0;
        for (int i = 0; i < gb->count; i++) {
            const rel_body_vector_t bv = craft_calcVectorToBody(sc, gb, craft_id, i);
            if (bv.distance < min_distance) {
                min_distance = bv.distance;
                closest_body_id = i;
            }
        }
        rel_body_vector_t nearest = craft_calcVectorToBody(sc, gb, craft_id, closest_body_id);

        snprintf(text_buffer, sizeof(text_buffer), "Nearest: %s", gb->names[closest_body_id]);
        SDL_WriteText(renderer, g_font_small, text_buffer, x_pos, y_offset, TEXT_COLOR);
        y_offset += line_height;

        snprintf(text_buffer, sizeof(text_buffer), "Distance: %.2f km", nearest.distance / 1000.0);
        SDL_WriteText(renderer, g_font_small, text_buffer, x_pos, y_offset, TEXT_COLOR);
        y_offset += line_height;

        snprintf(text_buffer, sizeof(text_buffer), "Rel. vel: %.2f m/s", nearest.rel_v);
        SDL_WriteText(renderer, g_font_small, text_buffer, x_pos, y_offset, TEXT_COLOR);
    }
}

// toggleable craft view window
void craft_RenderCraftView(SDL_Renderer* renderer, window_params_t* wp, body_properties_t* gb, spacecraft_properties_t* sc) {
    craft_DrawCraft(renderer, wp, sc, 0);
    craft_DrawRelativeArrows(renderer, wp, gb, sc, 0);
    craft_drawCraftPropertyArrows(renderer, wp, gb, sc, 0);
    craft_drawCraftInfo(renderer, wp, gb, sc, 0);
}
