#include "../gui/renderer.h"
#include "../globals.h"
#include "../sim/bodies.h"
#include "../sim/simulation.h"
#include "craft_view.h"
#include <math.h>
#include <stdio.h>

void readSimulationJSON(const char* FILENAME, body_properties_t* gb, spacecraft_properties_t* sc);

// initialize the window parameters
void init_window_params(window_params_t* wp) {
    // gets screen width and height parameters from computer
    const SDL_DisplayMode *mode = SDL_GetCurrentDisplayMode(SDL_GetPrimaryDisplay());

    wp->time_step = 1;
    // sets the default window size scaled based on the user's screen size
    wp->window_size_x = (float)mode->w * (2.0f/3.0f);
    wp->window_size_y = (float)mode->h * (2.0f/3.0f);
    wp->screen_origin_x = wp->window_size_x / 2;
    wp->screen_origin_y = wp->window_size_y / 2;
    wp->meters_per_pixel = 100000;
    wp->font_size = (float)wp->window_size_x / 50;
    wp->window_open = true;
    wp->sim_running = true;
    wp->data_logging_enabled = false;
    wp->sim_time = 0;

    wp->is_dragging = false;
    wp->drag_start_x = 0;
    wp->drag_start_y = 0;
    wp->drag_origin_x = 0;
    wp->drag_origin_y = 0;

    wp->main_view_shown = true;
    wp->craft_view_shown = false;
}

// draw a circle in SDL
void SDL_RenderFillCircle(SDL_Renderer* renderer, const float centerX, const float centerY, const float radius) {
    float x = radius;
    float y = 0;
    float radiusError = 1 - x;

    while (x >= y) {
        SDL_RenderLine(renderer, centerX - x, centerY + y, centerX + x, centerY + y);
        SDL_RenderLine(renderer, centerX - x, centerY - y, centerX + x, centerY - y);
        SDL_RenderLine(renderer, centerX - y, centerY + x, centerX + y, centerY + x);
        SDL_RenderLine(renderer, centerX - y, centerY - x, centerX + y, centerY - x);

        y++;
        if (radiusError < 0) {
            radiusError += 2 * y + 1;
        } else {
            x--;
            radiusError += 2 * (y - x + 1);
        }
    }
}

// helper function to write text to the screen
void SDL_WriteText(SDL_Renderer* renderer, TTF_Font* font, const char* text, const float x, const float y, const SDL_Color color) {
    if (!text || !font || !renderer) return;

    SDL_Surface* text_surface = TTF_RenderText_Blended(font, text, 0, color);
    if (!text_surface) return;

    SDL_Texture* text_texture = SDL_CreateTextureFromSurface(renderer, text_surface);
    if (!text_texture) {
        SDL_DestroySurface(text_surface);
        return;
    }

    SDL_FRect text_rect = {roundf(x), roundf(y), (float)text_surface->w, (float)text_surface->h};
    SDL_RenderTexture(renderer, text_texture, NULL, &text_rect);

    SDL_DestroyTexture(text_texture);
    SDL_DestroySurface(text_surface);
}

// draw graph-like axes with tick marks and labels
void drawScaleBar(SDL_Renderer* renderer, const window_params_t wp) {
    SDL_SetRenderDrawColor(renderer, 60, 70, 85, 255);

    float axis_thickness = 1.5f;
    float tick_length = 12;
    float tick_spacing_pixels = wp.window_size_x * 0.15f; // spacing between tick marks

    // calculate the distance each tick represents in meters
    float distance_per_tick = tick_spacing_pixels * (float)wp.meters_per_pixel;

    // determine magnitude and round to nice value
    float magnitude = powf(10.0f, floorf(log10f(distance_per_tick)));
    float nice_distance = roundf(distance_per_tick / magnitude) * magnitude;

    // recalculate pixel spacing for nice rounded values
    float nice_tick_spacing = nice_distance / (float)wp.meters_per_pixel;

    // draw X axis
    SDL_FRect x_axis = {0, wp.screen_origin_y - axis_thickness/2.0f, wp.window_size_x, axis_thickness};
    SDL_RenderFillRect(renderer, &x_axis);

    // draw Y axis
    SDL_FRect y_axis = {wp.screen_origin_x - axis_thickness/2.0f, 0, axis_thickness, wp.window_size_y};
    SDL_RenderFillRect(renderer, &y_axis);

    // draw tick marks and labels on X axis
    for (int i = -10; i <= 10; i++) {
        if (i == 0) continue; // skip center

        float tick_x = wp.screen_origin_x + (float)i * nice_tick_spacing;
        if (tick_x < 0 || tick_x > wp.window_size_x) continue;

        // draw tick mark
        SDL_FRect tick = {tick_x - 1, wp.screen_origin_y - tick_length/2.0f, 2, tick_length};
        SDL_RenderFillRect(renderer, &tick);

        // draw label
        double distance_value = (float)i * nice_distance;
        char label[50];
        if (fabs(distance_value) >= 1000000) {
            snprintf(label, sizeof(label), "%.0f km", distance_value / 1000.0);
        } else {
            snprintf(label, sizeof(label), "%.0f m", distance_value);
        }
        SDL_WriteText(renderer, g_font_small, label, tick_x - 20, wp.screen_origin_y + tick_length, (SDL_Color){140, 150, 165, 255});
    }

    // draw tick marks and labels on Y axis
    for (int i = -10; i <= 10; i++) {
        if (i == 0) continue; // skip center

        float tick_y = wp.screen_origin_y - (float)i * nice_tick_spacing;
        if (tick_y < 0 || tick_y > wp.window_size_y) continue;

        // draw tick mark
        SDL_FRect tick = {wp.screen_origin_x - tick_length/2.0f, tick_y - 1, tick_length, 2};
        SDL_RenderFillRect(renderer, &tick);

        // draw label
        float distance_value = (float)i * nice_distance;
        char label[50];
        if (fabsf(distance_value) >= 1000000) {
            snprintf(label, sizeof(label), "%.0f km", distance_value / 1000.0);
        } else {
            snprintf(label, sizeof(label), "%.0f m", distance_value);
        }
        SDL_WriteText(renderer, g_font_small, label, wp.screen_origin_x + tick_length, tick_y - 5, (SDL_Color){140, 150, 165, 255});
    }

    // draw origin label
    SDL_WriteText(renderer, g_font_small, "0", wp.screen_origin_x + 5, wp.screen_origin_y + 5, (SDL_Color){140, 150, 165, 255});
}

// thing that calculates changing sim speed
bool isMouseInRect(const int mouse_x, const int mouse_y, const int rect_x, const int rect_y, const int rect_w, const int rect_h) {
    return (mouse_x >= rect_x && mouse_x <= rect_x + rect_w &&
            mouse_y >= rect_y && mouse_y <= rect_y + rect_h);
}

void body_renderOrbitBodies(SDL_Renderer* renderer, body_properties_t* gb, window_params_t* wp) {
    for (int i = 0; i < gb->count; i++) {
        SDL_SetRenderDrawColor(renderer, 255, 255, 255, 255);
        // draw bodies
        SDL_RenderFillCircle(renderer, gb->pixel_coordinates_x[i],
                        gb->pixel_coordinates_y[i],
                        body_calculateVisualRadius(gb, i, *wp));
        SDL_WriteText(renderer, g_font_small, gb->names[i], gb->pixel_coordinates_x[i] - gb->pixel_radius[i], gb->pixel_coordinates_y[i] + gb->pixel_radius[i], TEXT_COLOR);

        // cache current position
        gb->cached_body_coords_x[i][gb->cache_counter] = gb->pixel_coordinates_x[i];
        gb->cached_body_coords_y[i][gb->cache_counter] = gb->pixel_coordinates_y[i];

        // if zooming, clear the path cache for this body
        if (wp->is_zooming || wp->is_dragging) {
            for (int n = 0; n < PATH_CACHE_LENGTH; n++) {
                gb->cached_body_coords_x[i][n] = gb->pixel_coordinates_x[i];
                gb->cached_body_coords_y[i][n] = gb->pixel_coordinates_y[i];
            }
        }

        // draw the paths of the bodies
        SDL_SetRenderDrawColor(renderer, ACCENT_COLOR.r, ACCENT_COLOR.g, ACCENT_COLOR.b, ACCENT_COLOR.a);
        for (int j = 0; j < PATH_CACHE_LENGTH; j++) {
            SDL_RenderPoint(renderer, gb->cached_body_coords_x[i][j], gb->cached_body_coords_y[i][j]);
        }
    }

    // reset flag after clearing paths for all bodies
    if (wp->is_zooming) wp->is_zooming = false;

    if (gb->cache_counter >= PATH_CACHE_LENGTH) gb->cache_counter = 0;
    gb->cache_counter++;
}

void craft_renderCrafts(SDL_Renderer* renderer, const spacecraft_properties_t* sc) {
    SDL_SetRenderDrawColor(renderer, 255, 100, 100, 255); // reddish color for spacecraft
    for (int i = 0; i < sc->count; i++) {
        if (sc->engine_on[i]) SDL_SetRenderDrawColor(renderer, 255, 190, 190, 255);
        // draw craft as a small filled circle
        SDL_RenderFillCircle(renderer, sc->pixel_coordinates_x[i], sc->pixel_coordinates_y[i], 3);
        // draw the name
        SDL_WriteText(renderer, g_font_small, sc->names[i], sc->pixel_coordinates_x[i] + 5, sc->pixel_coordinates_y[i] + 5, TEXT_COLOR);
    }
}

void renderTimeIndicators(SDL_Renderer* renderer, const window_params_t wp) {
    float right_margin = wp.window_size_x * 0.68f;
    float top_margin = wp.window_size_y * 0.015f;

    // show sim time in the top corner
    char time[32];
    if (wp.sim_time < 60) {
        snprintf(time, sizeof(time), "Time: %.f s", wp.sim_time); // sim time in seconds
    }
    else if (wp.sim_time > 60 && wp.sim_time < 3600) {
        snprintf(time, sizeof(time), "Time: %.f mins", wp.sim_time/60); // sim time in minutes
    }
    else if (wp.sim_time > 3600 && wp.sim_time < 86400) {
        snprintf(time, sizeof(time), "Time: %.1f hrs", wp.sim_time/3600); // sim time in hours
    }
    else {
        snprintf(time, sizeof(time), "Time: %.1f days", wp.sim_time/86400); // sim time in days
    }
    SDL_WriteText(renderer, g_font, time, right_margin, top_margin + wp.font_size * 1.2f, TEXT_COLOR);

    // show paused indication
    SDL_Color status_color = wp.sim_running ? (SDL_Color){120, 220, 140, 255} : (SDL_Color){240, 140, 120, 255};
    const char* status_text = wp.sim_running ? "Running" : "Paused";
    SDL_WriteText(renderer, g_font, status_text, right_margin, top_margin, status_color);

}

// creates a texture for a button with the given text and background color
SDL_Texture* createButtonTexture(SDL_Renderer* renderer, const button_t* button, const char* text, SDL_Color bg_color) {
    SDL_Texture* tex = SDL_CreateTexture(renderer, SDL_PIXELFORMAT_RGBA8888, SDL_TEXTUREACCESS_TARGET, (int)button->width, (int)button->height);
    if (!tex) return NULL;

    SDL_SetTextureBlendMode(tex, SDL_BLENDMODE_NONE);
    SDL_SetRenderTarget(renderer, tex);

    SDL_SetRenderDrawColor(renderer, 0, 0, 0, 0);
    SDL_RenderClear(renderer);

    // drop shadow
    SDL_SetRenderDrawColor(renderer, 0, 0, 0, 128);
    SDL_FRect shadow_rect = {1, 1, button->width - 2, button->height - 2};
    SDL_RenderFillRect(renderer, &shadow_rect);

    // main background
    SDL_SetRenderDrawColor(renderer, bg_color.r, bg_color.g, bg_color.b, bg_color.a);
    SDL_FRect bg_rect = {1, 1, button->width - 2, button->height - 2};
    SDL_RenderFillRect(renderer, &bg_rect);

    // border/shadow
    SDL_Color border_light = {200, 200, 200, 255};
    SDL_Color border_dark = {30, 30, 30, 255};

    if (button->is_hovered) {
        SDL_SetRenderDrawColor(renderer, border_dark.r, border_dark.g, border_dark.b, border_dark.a);
    }
    else {
        SDL_SetRenderDrawColor(renderer, border_light.r, border_light.g, border_light.b, border_light.a);
    }
    SDL_RenderLine(renderer, 0, 0, button->width - 2, 0);
    SDL_RenderLine(renderer, 0, 0, 0, button->height - 2);

    if (button->is_hovered) {
        SDL_SetRenderDrawColor(renderer, border_light.r, border_light.g, border_light.b, border_light.a);
    }
    else {
        SDL_SetRenderDrawColor(renderer, border_dark.r, border_dark.g, border_dark.b, border_dark.a);
    }
    SDL_RenderLine(renderer, 1, button->height - 2, button->width - 2, button->height - 2);
    SDL_RenderLine(renderer, button->width - 2, 1, button->width - 2, button->height - 2);

    // centered text
    if (g_font && text) {
        SDL_Surface* text_surface = TTF_RenderText_Solid(g_font, text, 0, TEXT_COLOR);
        if (text_surface) {
            float text_x = roundf((button->width - (float)text_surface->w) / 2.0f);
            float text_y = roundf((button->height - (float)text_surface->h) / 2.0f);

            SDL_Texture* text_tex = SDL_CreateTextureFromSurface(renderer, text_surface);
            if (text_tex) {
                SDL_FRect text_rect = {text_x, text_y, (float)text_surface->w, (float)text_surface->h};
                SDL_RenderTexture(renderer, text_tex, NULL, &text_rect);
                SDL_DestroyTexture(text_tex);
            }
            SDL_DestroySurface(text_surface);
        }
    }

    SDL_SetRenderTarget(renderer, NULL);
    return tex;
}



// renders a button using the default button texture
void renderButton(SDL_Renderer* renderer, const button_t* button) {
    // choose which texture to use based on hover state
    SDL_Texture* tex = button->is_hovered ? button->hover_texture : button->normal_texture;

    if (tex) {
        SDL_FRect dst = {button->x, button->y, button->width, button->height};
        SDL_RenderTexture(renderer, tex, NULL, &dst);
    }
}

// updates button textures (call when text changes or button is created)
void updateButtonTextures(SDL_Renderer* renderer, button_t* button, const char* text) {
    // destroy old textures if they exist
    if (button->normal_texture) {
        SDL_DestroyTexture(button->normal_texture);
        button->normal_texture = NULL;
    }
    if (button->hover_texture) {
        SDL_DestroyTexture(button->hover_texture);
        button->hover_texture = NULL;
    }

    // create new textures
    button->normal_texture = createButtonTexture(renderer, button, text, button->normal_color);
    button->hover_texture = createButtonTexture(renderer, button, text, button->hover_color);
    button->textures_valid = true;
}

// destroys button textures (for cleanup)
void destroyButtonTextures(button_t* button) {
    if (button->normal_texture) {
        SDL_DestroyTexture(button->normal_texture);
        button->normal_texture = NULL;
    }
    if (button->hover_texture) {
        SDL_DestroyTexture(button->hover_texture);
        button->hover_texture = NULL;
    }
    button->textures_valid = false;
}

// destroys all button textures in button_storage_t
void destroyAllButtonTextures(button_storage_t* buttons) {
    destroyButtonTextures(&buttons->sc_button);
    destroyButtonTextures(&buttons->csv_load_button);
    destroyButtonTextures(&buttons->craft_view_button);
    destroyButtonTextures(&buttons->show_stats_button);
}

// initializes all button properties based on window parameters
void initButtons(button_storage_t* buttons, const window_params_t wp) {
    float margin = wp.window_size_x * 0.015f;
    float button_spacing = wp.window_size_y * 0.01f;
    float button_height = wp.window_size_y * 0.045f;
    float side_button_width = wp.window_size_x * 0.14f;

    // speed control button (top left)
    buttons->sc_button = (button_t){
        .x = margin,
        .y = margin,
        .width = wp.window_size_x * 0.26f,
        .height = button_height,
        .is_hovered = false,
        .normal_color = BUTTON_COLOR,
        .hover_color = BUTTON_HOVER_COLOR,
        .normal_texture = NULL,
        .hover_texture = NULL,
        .textures_valid = false
    };

    // csv loading button (bottom right)
    buttons->csv_load_button = (button_t){
        .width = side_button_width,
        .height = button_height,
        .x = wp.window_size_x - side_button_width - margin,
        .y = wp.window_size_y - button_height - margin,
        .is_hovered = false,
        .normal_color = BUTTON_COLOR,
        .hover_color = BUTTON_HOVER_COLOR,
        .normal_texture = NULL,
        .hover_texture = NULL,
        .textures_valid = false
    };

    // craft view button (above csv button)
    buttons->craft_view_button = (button_t){
        .x = buttons->csv_load_button.x,
        .y = buttons->csv_load_button.y - button_height - button_spacing,
        .width = side_button_width,
        .height = button_height,
        .is_hovered = false,
        .normal_color = BUTTON_COLOR,
        .hover_color = BUTTON_HOVER_COLOR,
        .normal_texture = NULL,
        .hover_texture = NULL,
        .textures_valid = false
    };

    // show stats window button (above craft view button)
    buttons->show_stats_button = (button_t){
        .x = buttons->craft_view_button.x,
        .y = buttons->craft_view_button.y - button_height - button_spacing,
        .width = side_button_width,
        .height = button_height,
        .is_hovered = false,
        .normal_color = BUTTON_COLOR,
        .hover_color = BUTTON_HOVER_COLOR,
        .normal_texture = NULL,
        .hover_texture = NULL,
        .textures_valid = false
    };
}

// renders all the buttons on the screen, this function holds all button drawing logic
// remember: you must add this button to the event handling logic for it to work!!!
// remember: you must add each button to the button_storage_t struct!!
// remember: you must add this to initButtons!
// remember: you must add this to destroyAllButtonTextures
void renderUIButtons(SDL_Renderer* renderer, button_storage_t* buttons, const window_params_t* wp) {
    // speed control button - dynamic text, update every frame
    char speed_text[32];
    snprintf(speed_text, sizeof(speed_text), "Time Step: %.5f s", wp->time_step);
    updateButtonTextures(renderer, &buttons->sc_button, speed_text);
    renderButton(renderer, &buttons->sc_button);

    // buttons
    updateButtonTextures(renderer, &buttons->csv_load_button, "Load CSV");
    renderButton(renderer, &buttons->csv_load_button);

    updateButtonTextures(renderer, &buttons->craft_view_button, "Craft View");
    renderButton(renderer, &buttons->craft_view_button);

    updateButtonTextures(renderer, &buttons->show_stats_button, "Stats");
    renderButton(renderer, &buttons->show_stats_button);
}

// show error window
void displayError(const char* title, const char* message) {
    SDL_ShowSimpleMessageBox(
        SDL_MESSAGEBOX_INFORMATION,
        title,
        message,
        NULL
    );
}

// shows FPS
void showFPS(SDL_Renderer* renderer, const Uint64 frame_start_time, const Uint64 perf_freq, const window_params_t wp, const bool FPS_SHOWN) {
    const float target_frame_time = 1.0f / 60.0f; // 1/fps (in seconds)
    Uint64 frame_end = SDL_GetPerformanceCounter();
    float frame_time = (float)(frame_end - frame_start_time) / (float)perf_freq;

    // delay to limit FPS
    if (frame_time < target_frame_time) {
        SDL_Delay((Uint32)((target_frame_time - frame_time) * 1000.0f));
    }

    frame_end = SDL_GetPerformanceCounter();
    frame_time = (float)(frame_end - frame_start_time) / (float)perf_freq;

    // display the FPS
    if (FPS_SHOWN) {
        const float fps_value = 1.0f / frame_time;
        char fps[25];
        snprintf(fps, sizeof(fps), "%.1f FPS", fps_value);
        const SDL_Color fps_color = fps_value >= 55.0f ? (SDL_Color){120, 220, 140, 255} : (SDL_Color){240, 200, 120, 255};
        SDL_WriteText(renderer, g_font_small, fps, wp.window_size_x * 0.015f, wp.window_size_y - wp.font_size * 1.8f, fps_color);
    }
}

// the stats box that shows stats yay
void renderStatsBox(SDL_Renderer* renderer, body_properties_t* bodies, const spacecraft_properties_t* sc, const window_params_t wp, stats_window_t* stats_window) {
    const float margin_x    = (wp.window_size_x * 0.02f);
    const float top_y       = (wp.window_size_y * 0.1f);
    const float line_height = (wp.window_size_y * 0.025f);

    if (!bodies || bodies->count == 0) return;

    // update cache every several frames (whatever its set to)
    if (stats_window->frame_counter == 0) {
        // cache velocity text for all bodies
        stats_window->cached_body_count = bodies->count;
        for (int i = 0; i < bodies->count && i < 10; i++) {
            snprintf(stats_window->vel_text[i], sizeof(stats_window->vel_text[i]),
                     "Vel %s: %.2e m/s", bodies->names[i], bodies->vel[i]);
        }

        // calculate and cache total system energy
        double total_energy = calculateTotalSystemEnergy(bodies, sc);
        snprintf(stats_window->total_energy_text, sizeof(stats_window->total_energy_text),
                 "Total Energy: %.2e J", total_energy);

        // cache change in total energy
        snprintf(stats_window->delta_energy_text, sizeof(stats_window->delta_energy_text),
                 "E'(t): %.2e", total_energy - stats_window->previous_total_energy);
        stats_window->previous_total_energy = total_energy;

        // measure initial energy on first call or when explicitly reset
        if (!stats_window->measured_initial_energy) {
            stats_window->initial_total_energy = total_energy;
            stats_window->measured_initial_energy = true;
        }

        // cache error text and color
        double error;
        if (stats_window->initial_total_energy != 0) {
            error = fabs(stats_window->initial_total_energy - total_energy) / fabs(stats_window->initial_total_energy) * 100.0;
        } else {
            error = 0.0;
        }
        snprintf(stats_window->error_text, sizeof(stats_window->error_text), "Error: %.4f%%", error);
        stats_window->error_color = error < 0.1 ? (SDL_Color){120, 220, 140, 255} : (SDL_Color){240, 200, 120, 255};
    }

    // render cached text every frame
    float y = top_y;

    // render velocities from cache
    for (int i = 0; i < stats_window->cached_body_count && i < 10; i++) {
        SDL_WriteText(renderer, g_font_small, stats_window->vel_text[i], margin_x, y, (SDL_Color){180, 190, 210, 255});
        y += line_height;
    }

    // render total energy from cache
    SDL_WriteText(renderer, g_font_small, stats_window->total_energy_text, margin_x, y, ACCENT_COLOR);
    y += line_height;

    // render delta energy from cache
    SDL_WriteText(renderer, g_font_small, stats_window->delta_energy_text, margin_x, y, ACCENT_COLOR);
    y += line_height;

    // render error from cache
    SDL_WriteText(renderer, g_font_small, stats_window->error_text, margin_x, y, stats_window->error_color);

    // increment frame counter
    stats_window->frame_counter++;
    // update every n frames
    if (stats_window->frame_counter >= 15) {
        stats_window->frame_counter = 0;
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// EVENT HANDLING HELPER FUNCTIONS
////////////////////////////////////////////////////////////////////////////////////////////////////
// handles mouse motion events (dragging and button hover states)
static void handleMouseMotionEvent(const SDL_Event* event, window_params_t* wp, button_storage_t* buttons) {
    int mouse_x = (int)event->motion.x;
    int mouse_y = (int)event->motion.y;

    // viewport dragging
    if (wp->is_dragging) {
        int delta_x = mouse_x - (int)wp->drag_start_x;
        int delta_y = mouse_y - (int)wp->drag_start_y;
        wp->screen_origin_x = wp->drag_origin_x + (float)delta_x;
        wp->screen_origin_y = wp->drag_origin_y + (float)delta_y;
    }//

    // update hover state for speed control button
    buttons->sc_button.is_hovered = isMouseInRect(mouse_x, mouse_y,
        (int)buttons->sc_button.x, (int)buttons->sc_button.y,
        (int)buttons->sc_button.width, (int)buttons->sc_button.height);
    // update hover state for csv loading button
    buttons->csv_load_button.is_hovered = isMouseInRect(mouse_x, mouse_y,
        (int)buttons->csv_load_button.x, (int)buttons->csv_load_button.y,
        (int)buttons->csv_load_button.width, (int)buttons->csv_load_button.height);
    // update hover state for stats button
    buttons->show_stats_button.is_hovered = isMouseInRect(mouse_x, mouse_y,
        (int)buttons->show_stats_button.x, (int)buttons->show_stats_button.y,
        (int)buttons->show_stats_button.width, (int)buttons->show_stats_button.height);
    // update hover state for craft view button
    buttons->craft_view_button.is_hovered = isMouseInRect(mouse_x, mouse_y,
        (int)buttons->craft_view_button.x, (int)buttons->craft_view_button.y,
        (int)buttons->craft_view_button.width, (int)buttons->craft_view_button.height);
}

// handles mouse button down events (button clicks and drag start)
static void handleMouseButtonDownEvent(const SDL_Event* event, window_params_t* wp, body_properties_t* gb, spacecraft_properties_t* sc, const button_storage_t* buttons, stats_window_t* stats_window) {
    // check if right mouse button or middle mouse button (for dragging)
    if (event->button.button == SDL_BUTTON_RIGHT || event->button.button == SDL_BUTTON_MIDDLE) {
        wp->is_dragging = true;
        wp->drag_start_x = event->button.x;
        wp->drag_start_y = event->button.y;
        wp->drag_origin_x = wp->screen_origin_x;
        wp->drag_origin_y = wp->screen_origin_y;
    }
    // existing button click handling
    else if (buttons->csv_load_button.is_hovered) {
        // reads the JSON file associated with loading orbital bodies and spacecraft
        readSimulationJSON("simulation_data.json", gb, sc);
        wp->sim_time = 0;
        // reset energy measurement so stats recalculate with new bodies
        stats_window->measured_initial_energy = false;
    }
    else if(buttons->show_stats_button.is_hovered) {
        // toggle stats display in main window
        stats_window->is_shown = !stats_window->is_shown;
    }
    else if (buttons->craft_view_button.is_hovered) {
        // toggle zoomed in view of craft
        if (!wp->craft_view_shown) {
            wp->craft_view_shown = true;
        }
        else {
            wp->craft_view_shown = false;
        }
    }
}

// handles mouse button release events
static void handleMouseButtonUpEvent(const SDL_Event* event, window_params_t* wp) {
    if (event->button.button == SDL_BUTTON_RIGHT || event->button.button == SDL_BUTTON_MIDDLE) {
        wp->is_dragging = false;
    }
}

// handles mouse wheel events (zooming and speed control)
static void handleMouseWheelEvent(const SDL_Event* event, window_params_t* wp, const button_storage_t* buttons) {
    int mouse_x = (int)event->wheel.mouse_x;
    int mouse_y = (int)event->wheel.mouse_y;

    // if it's in the speed control button area
    if (isMouseInRect(mouse_x, mouse_y,
        (int)buttons->sc_button.x, (int)buttons->sc_button.y,
        (int)buttons->sc_button.width, (int)buttons->sc_button.height)) {

        if (event->wheel.y > 0) {
            wp->time_step *= 1.05;
        } else if (event->wheel.y < 0) {
            wp->time_step /= 1.05;
        }
    }
    else {
        if (event->wheel.y > 0) {
            wp->is_zooming = true;
            wp->meters_per_pixel *= 1.05;
        } else if (event->wheel.y < 0) {
            wp->is_zooming = true;
            wp->meters_per_pixel /= 1.05;
        }
    }
}

// handles keyboard events (pause/play, reset)
static void handleKeyboardEvent(const SDL_Event* event, window_params_t* wp) {
    if(event->key.key == SDLK_SPACE) {
        if (wp->sim_running == false) {
            wp->sim_running = true;
        }
        else if (wp->sim_running == true) {
            wp->sim_running = false;
        }
    }
    else if (event->key.key == SDLK_R) {
        wp->reset_sim = true;
    }
}

// handles window resize events
static void handleWindowResizeEvent(const SDL_Event* event, window_params_t* wp, button_storage_t* buttons) {
    // store old font size to check if it changed
    static float last_font_size = 0.0f;
    static float last_small_font_size = 0.0f;

    wp->window_size_x = (float)event->window.data1;
    wp->window_size_y = (float)event->window.data2;
    wp->screen_origin_x = wp->window_size_x / 2;
    wp->screen_origin_y = wp->window_size_y / 2;

    // calculate new font sizes with minimum size validation
    wp->font_size = wp->window_size_x / 75;
    if (wp->font_size < 8.0f) wp->font_size = 8.0f;  // minimum font size

    float small_font_size = (float)wp->window_size_x / 90;
    if (small_font_size < 6.0f) small_font_size = 6.0f;  // minimum small font size

    // only reload fonts if size actually changed (prevents rapid reload during resize)
    if (fabsf(wp->font_size - last_font_size) > 0.5f) {
        TTF_Font* new_font = TTF_OpenFont("font.ttf", wp->font_size);
        if (!new_font) {
            fprintf(stderr, "Failed to load font during resize (size: %.1f): %s\n", wp->font_size, SDL_GetError());
            // keep the old font if new one fails to load
        } else {
            // only close old font after successfully loading new one
            if (g_font) TTF_CloseFont(g_font);
            g_font = new_font;
            TTF_SetFontHinting(g_font, TTF_HINTING_NORMAL);
            last_font_size = wp->font_size;
        }
    }

    if (fabsf(small_font_size - last_small_font_size) > 0.5f) {
        TTF_Font* new_font_small = TTF_OpenFont("font.ttf", small_font_size);
        if (!new_font_small) {
            fprintf(stderr, "Failed to load small font during resize (size: %.1f): %s\n", small_font_size, SDL_GetError());
            // keep the old font if new one fails to load
        } else {
            // only close old font after successfully loading new one
            if (g_font_small) TTF_CloseFont(g_font_small);
            g_font_small = new_font_small;
            TTF_SetFontHinting(g_font_small, TTF_HINTING_NORMAL);
            last_small_font_size = small_font_size;
        }
    }

    // destroy old button textures before reinitializing
    destroyAllButtonTextures(buttons);

    // update button position and size
    initButtons(buttons, *wp);
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// MAIN EVENT CHECKING FUNCTION
////////////////////////////////////////////////////////////////////////////////////////////////////
// the event handling code... checks if events are happening for input and does a task based on that input
void runEventCheck(SDL_Event* event, window_params_t* wp, body_properties_t* gb, spacecraft_properties_t* sc, button_storage_t* buttons, stats_window_t* stats_window) {
    while (SDL_PollEvent(event)) {
        // check if x button is pressed to quit
        if (event->type == SDL_EVENT_QUIT) {
            wp->reset_sim = true;
            wp->window_open = false;
            wp->sim_running = false;
        }
        // check if mouse is moving to update hover state
        else if (event->type == SDL_EVENT_MOUSE_MOTION && event->window.windowID == wp->main_window_ID) {
            handleMouseMotionEvent(event, wp, buttons);
        }
        // check if mouse button is clicked
        else if (event->type == SDL_EVENT_MOUSE_BUTTON_DOWN && event->window.windowID == wp->main_window_ID) {
            handleMouseButtonDownEvent(event, wp, gb, sc, buttons, stats_window);
        }
        // check if mouse button is released
        else if (event->type == SDL_EVENT_MOUSE_BUTTON_UP && event->window.windowID == wp->main_window_ID) {
            handleMouseButtonUpEvent(event, wp);
        }
        // check if scroll
        else if (event->type == SDL_EVENT_MOUSE_WHEEL && event->window.windowID == wp->main_window_ID) {
            handleMouseWheelEvent(event, wp, buttons);
        }
        // check if keyboard key is pressed
        else if (event->type == SDL_EVENT_KEY_DOWN && event->window.windowID == wp->main_window_ID) {
            handleKeyboardEvent(event, wp);
        }
        // check if window is resized (only for main window)
        else if (event->type == SDL_EVENT_WINDOW_RESIZED &&
                 event->window.windowID == wp->main_window_ID) {
            handleWindowResizeEvent(event, wp, buttons);
        }
    }
}
