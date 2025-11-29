#include "sdl_elements.h"
#include "sim_calculations.h"
#include <math.h>
#include <stdio.h>
#include "config.h"

SDL_Color TEXT_COLOR = {240, 240, 245, 255};
SDL_Color BUTTON_COLOR = {45, 52, 70, 255};
SDL_Color BUTTON_HOVER_COLOR = {65, 75, 100, 255};
SDL_Color ACCENT_COLOR = {80, 150, 220, 255};

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
    wp->sim_time = 0;

    wp->is_dragging = false;
    wp->drag_start_x = 0;
    wp->drag_start_y = 0;
    wp->drag_origin_x = 0;
    wp->drag_origin_y = 0;
}

// draw a circle in SDL
void SDL_RenderFillCircle(SDL_Renderer* renderer, float centerX, float centerY, float radius) {
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
void SDL_WriteText(SDL_Renderer* renderer, TTF_Font* font, const char* text, float x, float y, SDL_Color color) {
    if (!text || !font || !renderer) return;
    
    SDL_Surface* text_surface = TTF_RenderText_Blended(font, text, 0, color);
    if (!text_surface) return;
    
    SDL_Texture* text_texture = SDL_CreateTextureFromSurface(renderer, text_surface);
    if (!text_texture) {
        SDL_DestroySurface(text_surface);
        return;
    }
    
    SDL_FRect text_rect = {x, y, (float)text_surface->w, (float)text_surface->h};
    SDL_RenderTexture(renderer, text_texture, NULL, &text_rect);
    
    SDL_DestroyTexture(text_texture);
    SDL_DestroySurface(text_surface);
}

// draw graph-like axes with tick marks and labels
void drawScaleBar(SDL_Renderer* renderer, window_params_t wp) {
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
bool isMouseInRect(int mouse_x, int mouse_y, int rect_x, int rect_y, int rect_w, int rect_h) {
    return (mouse_x >= rect_x && mouse_x <= rect_x + rect_w &&
            mouse_y >= rect_y && mouse_y <= rect_y + rect_h);
}

void body_renderOrbitBodies(SDL_Renderer* renderer, body_properties_t* gb, window_params_t wp) {
    SDL_SetRenderDrawColor(renderer, 255, 255, 255, 255);
    for (int i = 0; i < gb->count; i++) {
        // draw bodies
        SDL_RenderFillCircle(renderer, gb->pixel_coordinates_x[i],
                        gb->pixel_coordinates_y[i],
                        body_calculateVisualRadius(gb, i, wp));
        SDL_WriteText(renderer, g_font_small, gb->names[i], gb->pixel_coordinates_x[i] - gb->pixel_radius[i], gb->pixel_coordinates_y[i] + gb->pixel_radius[i], TEXT_COLOR);
    }
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

void renderTimeIndicators(SDL_Renderer* renderer, window_params_t wp) {
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

// generic button renderer
void renderButton(SDL_Renderer* renderer, const button_t* button, const char* text) {
    // set background color based on hover state
    SDL_Color bg_color = button->is_hovered ? button->hover_color : button->normal_color;

    // draw button background
    SDL_FRect bg_rect = {
        (float)button->x,
        (float)button->y,
        (float)button->width,
        (float)button->height
    };
    SDL_SetRenderDrawColor(renderer, bg_color.r, bg_color.g, bg_color.b, bg_color.a);
    SDL_RenderFillRect(renderer, &bg_rect);

    // center text in button
    if (g_font && text) {
        SDL_Surface* text_surface = TTF_RenderText_Blended(g_font, text, 0, TEXT_COLOR);
        if (text_surface) {
            float text_x = button->x + (button->width - (float)text_surface->w) / 2.0f;
            float text_y = button->y + (button->height - (float)text_surface->h) / 2.0f;

            SDL_Texture* text_texture = SDL_CreateTextureFromSurface(renderer, text_surface);
            if (text_texture) {
                SDL_FRect text_rect = {text_x, text_y, (float)text_surface->w, (float)text_surface->h};
                SDL_RenderTexture(renderer, text_texture, NULL, &text_rect);
                SDL_DestroyTexture(text_texture);
            }
            SDL_DestroySurface(text_surface);
        }
    }
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
        .width = wp.window_size_x * 0.22f,
        .height = button_height,
        .is_hovered = false,
        .normal_color = BUTTON_COLOR,
        .hover_color = BUTTON_HOVER_COLOR
    };

    // csv loading button (bottom right)
    buttons->csv_load_button = (button_t){
        .width = side_button_width,
        .height = button_height,
        .x = wp.window_size_x - side_button_width - margin,
        .y = wp.window_size_y - button_height - margin,
        .is_hovered = false,
        .normal_color = BUTTON_COLOR,
        .hover_color = BUTTON_HOVER_COLOR
    };

    // show stats window button (above csv button)
    buttons->show_stats_button = (button_t){
        .x = buttons->csv_load_button.x,
        .y = buttons->csv_load_button.y - button_height - button_spacing,
        .width = side_button_width,
        .height = button_height,
        .is_hovered = false,
        .normal_color = BUTTON_COLOR,
        .hover_color = BUTTON_HOVER_COLOR
    };
}

// renders all the buttons on the screen, this function holds all button drawing logic
// remember: you must add this button to the event handling logic for it to work!!!
// remember: you must add each button to the button_storage_t struct!
void renderUIButtons(SDL_Renderer* renderer, const button_storage_t* buttons, const window_params_t* wp) {
    // speed control button
    char speed_text[32];
    snprintf(speed_text, sizeof(speed_text), "Time Step: %.2f s", wp->time_step);
    renderButton(renderer, &buttons->sc_button, speed_text);

    // csv loading button
    renderButton(renderer, &buttons->csv_load_button, "Load CSV");

    // show stats window button
    renderButton(renderer, &buttons->show_stats_button, "Stats");
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

void eror() {
    SDL_ShowSimpleMessageBox(
        SDL_MESSAGEBOX_INFORMATION,
        "hello",
        "you arem fuck",
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
void renderStatsBox(SDL_Renderer* renderer, body_properties_t* bodies, const spacecraft_properties_t* sc, window_params_t wp, stats_window_t* stats_window) {
    const float margin_x    = (wp.window_size_x * 0.02f);
    const float top_y       = (wp.window_size_y * 0.1f);
    const float line_height = (wp.window_size_y * 0.025f);

    if (!bodies || bodies->count == 0) return;

    float y = top_y;

    // velocities
    for (int i = 0; i < bodies->count; i++) {
        char vel_text[64];
        snprintf(vel_text, sizeof(vel_text), "Vel %s: %.2e m/s", bodies->names[i], bodies->vel[i]);
        SDL_WriteText(renderer, g_font_small, vel_text, margin_x, y, (SDL_Color){180, 190, 210, 255});
        y += line_height;
    }

    // total system energy
    double total_energy = calculateTotalSystemEnergy(bodies, sc);
    char total_energy_text[64];
    snprintf(total_energy_text, sizeof(total_energy_text), "Total Energy: %.2e J", total_energy);
    SDL_WriteText(renderer, g_font_small, total_energy_text, margin_x, y, ACCENT_COLOR);
    y += line_height;

    // total error - measure initial energy on first call or when explicitly reset
    if (!stats_window->measured_initial_energy) {
        stats_window->initial_total_energy = total_energy;
        stats_window->measured_initial_energy = true;
    }
    double error = (stats_window->initial_total_energy != 0) ? fabs(stats_window->initial_total_energy - total_energy) / fabs(stats_window->initial_total_energy) * 100.0 : 0.0;
    char error_text[64];
    snprintf(error_text, sizeof(error_text), "Error: %.4f%%", error);
    SDL_Color error_color = error < 0.1 ? (SDL_Color){120, 220, 140, 255} : (SDL_Color){240, 200, 120, 255};
    SDL_WriteText(renderer, g_font_small, error_text, margin_x, y, error_color);
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
    }

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
        // reads the JSON file associated with loading orbital bodies
        readBodyJSON("planet_data.json", gb);
        readSpacecraftJSON("spacecraft_data.json", sc);
        wp->sim_time = 0;
        // reset energy measurement so stats recalculate with new bodies
        stats_window->measured_initial_energy = false;
    }
    else if(buttons->show_stats_button.is_hovered) {
        // toggle stats display in main window
        stats_window->is_shown = !stats_window->is_shown;
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
            wp->meters_per_pixel *= 1.05;
        } else if (event->wheel.y < 0) {
            wp->meters_per_pixel /= 1.05;
        }
    }
}

// handles keyboard events (pause/play, reset)
static void handleKeyboardEvent(const SDL_Event* event, window_params_t* wp, body_properties_t* gb, spacecraft_properties_t* sc) {
    if(event->key.key == SDLK_SPACE) {
        if (wp->sim_running == false) {
            wp->sim_running = true;
        }
        else if (wp->sim_running == true) {
            wp->sim_running = false;
        }
    }
    else if (event->key.key == SDLK_R) {
        resetSim(&wp->sim_time, gb, sc);
    }
}

// handles window resize events
static void handleWindowResizeEvent(const SDL_Event* event, window_params_t* wp, button_storage_t* buttons) {
    wp->window_size_x = (float)event->window.data1;
    wp->window_size_y = (float)event->window.data2;
    wp->screen_origin_x = wp->window_size_x / 2;
    wp->screen_origin_y = wp->window_size_y / 2;

    // update font size proportionally
    wp->font_size = wp->window_size_x / 75;
    if (g_font) TTF_CloseFont(g_font);
    g_font = TTF_OpenFont("CascadiaCode.ttf", wp->font_size);

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
            handleKeyboardEvent(event, wp, gb, sc);
        }
        // check if window is resized (only for main window)
        else if (event->type == SDL_EVENT_WINDOW_RESIZED &&
                 event->window.windowID == wp->main_window_ID) {
            handleWindowResizeEvent(event, wp, buttons);
        }
    }
}