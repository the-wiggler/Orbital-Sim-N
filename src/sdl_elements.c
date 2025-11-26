#include "sdl_elements.h"
#include "sim_calculations.h"
#include <errno.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "config.h"

SDL_Color TEXT_COLOR = {255, 255, 255, 255};
SDL_Color BUTTON_COLOR = {50, 50, 50, 255};
SDL_Color BUTTON_HOVER_COLOR = {0, 0, 0, 255};

// initialize the window parameters
void init_window_params(window_params_t* wp) {
    wp->time_step = 1;
    wp->window_size_x = 1000;
    wp->window_size_y = 1000;
    wp->screen_origin_x = wp->window_size_x / 2;
    wp->screen_origin_y = wp->window_size_y / 2;
    wp->meters_per_pixel = 100000;
    wp->font_size = wp->window_size_x / 50;
    wp->window_open = true;
    wp->sim_running = true;
    wp->sim_time = 0;

    wp->is_dragging = false;
    wp->drag_start_x = 0;
    wp->drag_start_y = 0;
    wp->drag_origin_x = 0;
    wp->drag_origin_y = 0;
}

// initialize the text input dialogue stuff
void init_text_dialog(text_input_dialog_t* dialog) {
    dialog->active = false;
    dialog->state = INPUT_NONE;
}

// draw a circle in SDL
void SDL_RenderFillCircle(SDL_Renderer* renderer, int centerX, int centerY, int radius) {
    int x = radius;
    int y = 0;
    int radiusError = 1 - x;

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

// helper functiont to write text to the screen
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
    SDL_SetRenderDrawColor(renderer, 30, 30, 30, 255);

    int axis_thickness = 1;
    int tick_length = 10;
    int tick_spacing_pixels = wp.window_size_x * 0.15; // spacing between tick marks

    // calculate the distance each tick represents in meters
    double distance_per_tick = tick_spacing_pixels * wp.meters_per_pixel;

    // determine magnitude and round to nice value
    double magnitude = pow(10.0, floor(log10(distance_per_tick)));
    double nice_distance = round(distance_per_tick / magnitude) * magnitude;

    // recalculate pixel spacing for nice rounded values
    int nice_tick_spacing = (int)(nice_distance / wp.meters_per_pixel);

    // draw X axis
    SDL_FRect x_axis = {0, wp.screen_origin_y - axis_thickness/2.0f, wp.window_size_x, axis_thickness};
    SDL_RenderFillRect(renderer, &x_axis);

    // draw Y axis
    SDL_FRect y_axis = {wp.screen_origin_x - axis_thickness/2.0f, 0, axis_thickness, wp.window_size_y};
    SDL_RenderFillRect(renderer, &y_axis);

    // draw tick marks and labels on X axis
    for (int i = -10; i <= 10; i++) {
        if (i == 0) continue; // skip center

        int tick_x = wp.screen_origin_x + i * nice_tick_spacing;
        if (tick_x < 0 || tick_x > wp.window_size_x) continue;

        // draw tick mark
        SDL_FRect tick = {tick_x - 1, wp.screen_origin_y - tick_length/2.0f, 2, tick_length};
        SDL_RenderFillRect(renderer, &tick);

        // draw label
        double distance_value = i * nice_distance;
        char label[50];
        if (fabs(distance_value) >= 1000000) {
            sprintf(label, "%.0f km", distance_value / 1000.0);
        } else {
            sprintf(label, "%.0f m", distance_value);
        }
        SDL_WriteText(renderer, g_font_small, label, tick_x - 20, wp.screen_origin_y + tick_length, (SDL_Color){90, 90, 90, 255});
    }

    // draw tick marks and labels on Y axis
    for (int i = -10; i <= 10; i++) {
        if (i == 0) continue; // skip center

        int tick_y = wp.screen_origin_y - i * nice_tick_spacing;
        if (tick_y < 0 || tick_y > wp.window_size_y) continue;

        // draw tick mark
        SDL_FRect tick = {wp.screen_origin_x - tick_length/2.0f, tick_y - 1, tick_length, 2};
        SDL_RenderFillRect(renderer, &tick);

        // draw label
        double distance_value = i * nice_distance;
        char label[50];
        if (fabs(distance_value) >= 1000000) {
            sprintf(label, "%.0f km", distance_value / 1000.0);
        } else {
            sprintf(label, "%.0f m", distance_value);
        }
        SDL_WriteText(renderer, g_font_small, label, wp.screen_origin_x + tick_length, tick_y - 5, (SDL_Color){90, 90, 90, 255});
    }

    // draw origin label
    SDL_WriteText(renderer, g_font_small, "0", wp.screen_origin_x + 5, wp.screen_origin_y + 5, (SDL_Color){90, 90, 90, 255});
}

// thing that calculates changing sim speed
bool isMouseInRect(int mouse_x, int mouse_y, int rect_x, int rect_y, int rect_w, int rect_h) {
    return (mouse_x >= rect_x && mouse_x <= rect_x + rect_w &&
            mouse_y >= rect_y && mouse_y <= rect_y + rect_h);
}

void body_renderOrbitBodies(SDL_Renderer* renderer, body_properties_t* gb, int num_bodies, window_params_t wp) {
    SDL_SetRenderDrawColor(renderer, 255, 255, 255, 255);
    for (int i = 0; i < num_bodies; i++) {
        // draw bodies
        SDL_RenderFillCircle(renderer, gb[i].pixel_coordinates_x,
                        gb[i].pixel_coordinates_y, 
                        body_calculateVisualRadius(&gb[i], wp));
        SDL_WriteText(renderer, g_font_small, gb[i].name, gb[i].pixel_coordinates_x - gb[i].pixel_radius, gb[i].pixel_coordinates_y + gb[i].pixel_radius, TEXT_COLOR);
    }
}

void craft_renderCrafts(SDL_Renderer* renderer, spacecraft_properties_t* sc, int num_craft, window_params_t wp) {
    SDL_SetRenderDrawColor(renderer, 255, 100, 100, 255); // reddish color for spacecraft
    for (int i = 0; i < num_craft; i++) {
        if (sc[i].engine_on) SDL_SetRenderDrawColor(renderer, 255, 190, 190, 255);
        // draw craft as a small filled circle
        SDL_RenderFillCircle(renderer, sc[i].pixel_coordinates_x, sc[i].pixel_coordinates_y, 3);
        // draw the name
        SDL_WriteText(renderer, g_font_small, sc[i].name, sc[i].pixel_coordinates_x + 5, sc[i].pixel_coordinates_y + 5, TEXT_COLOR);
    }
}

void renderTimeIndicators(SDL_Renderer* renderer, window_params_t wp) {
    // show sim time in the top corner
    char time[32];
    if (wp.sim_time < 60) {
        snprintf(time, sizeof(time), "Sim time: %.f s", wp.sim_time); // sim time in seconds
    }
    else if (wp.sim_time > 60 && wp.sim_time < 3600) {
        snprintf(time, sizeof(time), "Sim time: %.f mins", wp.sim_time/60); // sim time in minutes
    }
    else if (wp.sim_time > 3600 && wp.sim_time < 86400) {
        snprintf(time, sizeof(time), "Sim time: %.1f hrs", wp.sim_time/3600); // sim time in hours
    }
    else {
        snprintf(time, sizeof(time), "Sim time: %.1f days", wp.sim_time/86400); // sim time in days
    }
    SDL_WriteText(renderer, g_font, time, wp.window_size_x * 0.75, wp.window_size_y * 0.04, TEXT_COLOR);

    // show paused indication
    if (wp.sim_running) SDL_WriteText(renderer, g_font, "Sim Running...", wp.window_size_x * 0.75, wp.window_size_y * 0.015, TEXT_COLOR);
    else SDL_WriteText(renderer, g_font, "Sim Paused", wp.window_size_x * 0.75, wp.window_size_y * 0.015, TEXT_COLOR);

}

// generic button renderer
void renderButton(SDL_Renderer* renderer, button_t* button, const char* text, window_params_t wp) {
    // set background color based on hover state
    SDL_Color bg_color = button->is_hovered ? button->hover_color : button->normal_color;
    SDL_SetRenderDrawColor(renderer, bg_color.r, bg_color.g, bg_color.b, bg_color.a);

    // draw button background
    SDL_FRect bg_rect = {
        (float)button->x,
        (float)button->y,
        (float)button->width,
        (float)button->height
    };
    SDL_RenderFillRect(renderer, &bg_rect);

    // center text in button
    int padding_x = wp.window_size_x * 0.01;
    int padding_y = wp.window_size_y * 0.01;
    SDL_WriteText(renderer, g_font, text, button->x + padding_x, button->y + padding_y, TEXT_COLOR);
}

// initializes all button properties based on window parameters
void initButtons(button_storage_t* buttons, window_params_t wp) {
    // speed control button
    buttons->sc_button = (button_t){
        .x = wp.window_size_x * 0.01,
        .y = wp.window_size_y * 0.007,
        .width = wp.window_size_x * 0.25,
        .height = wp.window_size_y * 0.04,
        .is_hovered = false,
        .normal_color = BUTTON_COLOR,
        .hover_color = BUTTON_HOVER_COLOR
    };

    // csv loading button
    buttons->csv_load_button = (button_t){
        .width = wp.window_size_x * 0.12,
        .height = wp.window_size_y * 0.04,
        .x = wp.window_size_x - wp.window_size_x * 0.12 - 10,
        .y = wp.window_size_y - wp.window_size_y * 0.04 - 10,
        .is_hovered = false,
        .normal_color = BUTTON_COLOR,
        .hover_color = BUTTON_HOVER_COLOR
    };

    // add body button
    buttons->add_body_button = (button_t){
        .x = buttons->csv_load_button.x,
        .y = buttons->csv_load_button.y - buttons->csv_load_button.height,
        .width = buttons->csv_load_button.width,
        .height = buttons->csv_load_button.height,
        .is_hovered = false,
        .normal_color = BUTTON_COLOR,
        .hover_color = BUTTON_HOVER_COLOR
    };

    // show stats window button
    buttons->show_stats_button = (button_t){
        .x = buttons->add_body_button.x,
        .y = buttons->add_body_button.y - buttons->add_body_button.height,
        .width = buttons->add_body_button.width,
        .height = buttons->add_body_button.height,
        .is_hovered = false,
        .normal_color = BUTTON_COLOR,
        .hover_color = BUTTON_HOVER_COLOR
    };
}

// renders all of the buttons on the screen, this function holds all button drawing logic
// remember: you must add this button to the event handling logic for it to work!!!
// remember: you must add each button to the button_storage_t struct!
void renderUIButtons(SDL_Renderer* renderer, button_storage_t* buttons, window_params_t* wp) {
    // speed control button
    char speed_text[32];
    snprintf(speed_text, sizeof(speed_text), "Speed: %.2f s/frame", wp->time_step);
    renderButton(renderer, &buttons->sc_button, speed_text, *wp);

    // csv loading button
    renderButton(renderer, &buttons->csv_load_button, "Load CSV", *wp);

    // add orbital body button
    renderButton(renderer, &buttons->add_body_button, "Add Body", *wp);

    // show stats window button
    renderButton(renderer, &buttons->show_stats_button, "Stats", *wp);
}

// renders the text input dialog box
void renderBodyTextInputDialog(SDL_Renderer* renderer, text_input_dialog_t* dialog, window_params_t wp) {
    if (!dialog->active) return;
    // background overlay
    SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
    SDL_FRect overlay = {0, 0, (float)wp.window_size_x, (float)wp.window_size_y};
    SDL_RenderFillRect(renderer, &overlay);
    // dialog box dimensions
    int dialog_width = wp.window_size_x * 0.5;
    int dialog_height = wp.window_size_y * 0.3;
    int dialog_x = (wp.window_size_x - dialog_width) / 2;
    int dialog_y = (wp.window_size_y - dialog_height) / 2;
    // draw dialog background
    SDL_SetRenderDrawColor(renderer, 40, 40, 40, 255);
    SDL_FRect dialog_bg = {(float)dialog_x, (float)dialog_y, (float)dialog_width, (float)dialog_height};
    SDL_RenderFillRect(renderer, &dialog_bg);
    // draw dialog border
    SDL_SetRenderDrawColor(renderer, 255, 255, 255, 255);
    SDL_RenderRect(renderer, &dialog_bg);

    // determine prompt text based on current state
    const char* prompt_text;
    switch (dialog->state) {
        case INPUT_NAME:   prompt_text = "Enter name:"; break;
        case INPUT_MASS:   prompt_text = "Enter mass (kg):"; break;
        case INPUT_X_POS:  prompt_text = "Enter X position (m):"; break;
        case INPUT_Y_POS:  prompt_text = "Enter Y position (m):"; break;
        case INPUT_X_VEL:  prompt_text = "Enter X velocity (m/s):"; break;
        case INPUT_Y_VEL:  prompt_text = "Enter Y velocity (m/s):"; break;
        default:           prompt_text = ""; break;
    }

    // render prompt text
    int text_x = dialog_x + dialog_width * 0.1;
    int text_y = dialog_y + dialog_height * 0.2;
    SDL_WriteText(renderer, g_font, prompt_text, text_x, text_y, (SDL_Color){255, 255, 255, 255});

    // render current input with cursor
    char display_buffer[300];
    snprintf(display_buffer, sizeof(display_buffer), "%s_", dialog->input_buffer);
    SDL_WriteText(renderer, g_font, display_buffer, text_x, text_y + wp.font_size * 1.5, (SDL_Color){200, 100, 200, 255});

    // render instruction text
    const char* instruction = "Press Enter to continue | ESC to cancel";
    SDL_WriteText(renderer, g_font, instruction, text_x, text_y + wp.font_size * 3.5, (SDL_Color){180, 180, 180, 255});
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

// helper function to validate if a string is a valid number
bool isValidNumber(const char* str, double* out_value) {
    if (!str || strlen(str) == 0) {
        return false;
    }

    char* endptr;
    errno = 0;
    double value = strtod(str, &endptr);

    // check if conversion failed
    if (endptr == str || *endptr != '\0' || errno == ERANGE) {
        return false;
    }

    // check for NaN or infinity
    if (isnan(value) || isinf(value)) {
        return false;
    }

    if (out_value) {
        *out_value = value;
    }

    return true;
}

// shows FPS
void showFPS(SDL_Renderer* renderer, Uint64 frame_start_timem, Uint64 perf_freq, window_params_t wp) {
    // end frame and calculate FPS
    Uint64 frame_end = SDL_GetPerformanceCounter();
    double dt = (double)(frame_end - frame_start_timem) / (double)perf_freq;
    char fps[25];
    snprintf(fps, sizeof(fps), "%.1f FPS", 1.0 / dt);
    SDL_WriteText(renderer, g_font, fps, wp.window_size_x * 0.01, wp.window_size_y - 0.03 * wp.window_size_y, TEXT_COLOR);
}

// the stats box that shows stats yay
void renderStatsBox(SDL_Renderer* renderer, body_properties_t* bodies, int num_bodies, spacecraft_properties_t* sc, int num_craft, window_params_t wp, stats_window_t* stats_window) {
    int margin_x    = (int)(wp.window_size_x * 0.02f);
    int top_y       = (int)(wp.window_size_y * 0.07f);
    int line_height = (int)(wp.window_size_y * 0.02f);

    if (!bodies || num_bodies == 0) return;

    int y = top_y;

    // velocities
    for (int i = 0; i < num_bodies; i++) {
        char vel_text[64];
        snprintf(vel_text, sizeof(vel_text), "Vel of %s: %.2e m/s", bodies[i].name, bodies[i].vel);
        SDL_WriteText(renderer, g_font, vel_text, margin_x, y, TEXT_COLOR);
        y += line_height;
    }

    // kinetic energies
    for (int i = 0; i < num_bodies; i++) {
        char ke_text[64];
        snprintf(ke_text, sizeof(ke_text), "KE of %s: %.2e J", bodies[i].name, bodies[i].kinetic_energy);
        SDL_WriteText(renderer, g_font, ke_text, margin_x, y, TEXT_COLOR);
        y += line_height;
    }

    // total system energy
    double total_energy = calculateTotalSystemEnergy(bodies, sc, num_bodies, num_craft);
    char total_energy_text[64];
    snprintf(total_energy_text, sizeof(total_energy_text), "Total System Energy: %.2e J", total_energy);
    SDL_WriteText(renderer, g_font, total_energy_text, margin_x, y, TEXT_COLOR);
    y += line_height;

    // total error - measure initial energy on first call or when explicitly reset
    if (!stats_window->measured_initial_energy) {
        stats_window->initial_total_energy = total_energy;
        stats_window->measured_initial_energy = true;
    }
    double error = (stats_window->initial_total_energy != 0) ? fabs(stats_window->initial_total_energy - total_energy) / fabs(stats_window->initial_total_energy) * 100.0 : 0.0;
    char error_text[64];
    snprintf(error_text, sizeof(error_text), "Energy Error: %.4f%%", error);
    SDL_WriteText(renderer, g_font, error_text, margin_x, y, TEXT_COLOR);
}


////////////////////////////////////////////////////////////////////////////////////////////////////
// EVENT HANDLING HELPER FUNCTIONS
////////////////////////////////////////////////////////////////////////////////////////////////////
// handles text input dialog events
static bool handleTextInputDialogEvent(SDL_Event* event, window_params_t* wp, text_input_dialog_t* dialog, body_properties_t** gb, int* num_bodies) {
    if (event->type == SDL_EVENT_TEXT_INPUT && event->window.windowID == wp->main_window_ID) {
        // append character to input buffer
        size_t len = strlen(dialog->input_buffer);
        if (len < sizeof(dialog->input_buffer) - 1) {
            strncat(dialog->input_buffer, event->text.text, sizeof(dialog->input_buffer) - len - 1);
        }
        return true;
    }
    else if (event->type == SDL_EVENT_KEY_DOWN) {
        if (event->key.key == SDLK_BACKSPACE && strlen(dialog->input_buffer) > 0) {
            // remove last character
            dialog->input_buffer[strlen(dialog->input_buffer) - 1] = '\0';
            return true;
        }
        else if ((event->key.key == SDLK_RETURN || event->key.key == SDLK_KP_ENTER) && event->window.windowID == wp->main_window_ID) {
            // process current input and move to next state
            switch (dialog->state) {
                case INPUT_NAME:
                    if (strlen(dialog->input_buffer) == 0) {
                        displayError("Invalid Input", "Name cannot be empty");
                        break;
                    }
                    strncpy(dialog->name, dialog->input_buffer, sizeof(dialog->name) - 1);
                    dialog->state = INPUT_MASS;
                    break;
                case INPUT_MASS: {
                    double mass;
                    if (!isValidNumber(dialog->input_buffer, &mass)) {
                        displayError("Invalid Input", "Mass must be a valid number");
                        break;
                    }
                    if (mass <= 0) {
                        displayError("Invalid Input", "Mass must be positive (greater than 0)");
                        break;
                    }
                    dialog->mass = mass;
                    dialog->state = INPUT_X_POS;
                    break;
                }
                case INPUT_X_POS: {
                    double x_pos;
                    if (!isValidNumber(dialog->input_buffer, &x_pos)) {
                        displayError("Invalid Input", "X position must be a valid number");
                        break;
                    }
                    dialog->x_pos = x_pos;
                    dialog->state = INPUT_Y_POS;
                    break;
                }
                case INPUT_Y_POS: {
                    double y_pos;
                    if (!isValidNumber(dialog->input_buffer, &y_pos)) {
                        displayError("Invalid Input", "Y position must be a valid number");
                        break;
                    }
                    dialog->y_pos = y_pos;
                    dialog->state = INPUT_X_VEL;
                    break;
                }
                case INPUT_X_VEL: {
                    double x_vel;
                    if (!isValidNumber(dialog->input_buffer, &x_vel)) {
                        displayError("Invalid Input", "X velocity must be a valid number");
                        break;
                    }
                    dialog->x_vel = x_vel;
                    dialog->state = INPUT_Y_VEL;
                    break;
                }
                case INPUT_Y_VEL: {
                    double y_vel;
                    if (!isValidNumber(dialog->input_buffer, &y_vel)) {
                        displayError("Invalid Input", "Y velocity must be a valid number");
                        break;
                    }
                    dialog->y_vel = y_vel;
                    // all parameters collected, add the body
                    body_addOrbitalBody(gb, num_bodies, dialog->name, dialog->mass,
                             dialog->x_pos, dialog->y_pos, dialog->x_vel, dialog->y_vel);
                    // reset dialog
                    dialog->active = false;
                    dialog->state = INPUT_NONE;
                    SDL_StopTextInput(SDL_GetKeyboardFocus());
                    break;
                }
                default:
                    break;
            }
            // clear input buffer for next field
            dialog->input_buffer[0] = '\0';
            return true;
        }
        else if (event->key.key == SDLK_ESCAPE) {
            // cancel dialog
            dialog->active = false;
            dialog->state = INPUT_NONE;
            dialog->input_buffer[0] = '\0';
            SDL_StopTextInput(SDL_GetKeyboardFocus());
            return true;
        }
    }
    return false;
}

// handles mouse motion events (dragging and button hover states)
static void handleMouseMotionEvent(SDL_Event* event, window_params_t* wp, button_storage_t* buttons) {
    int mouse_x = (int)event->motion.x;
    int mouse_y = (int)event->motion.y;

    // viewport dragging
    if (wp->is_dragging) {
        int delta_x = mouse_x - wp->drag_start_x;
        int delta_y = mouse_y - wp->drag_start_y;
        wp->screen_origin_x = wp->drag_origin_x + delta_x;
        wp->screen_origin_y = wp->drag_origin_y + delta_y;
    }

    // update hover state for speed control button
    buttons->sc_button.is_hovered = isMouseInRect(mouse_x, mouse_y,
        buttons->sc_button.x, buttons->sc_button.y,
        buttons->sc_button.width, buttons->sc_button.height);
    // update hover state for csv loading button
    buttons->csv_load_button.is_hovered = isMouseInRect(mouse_x, mouse_y,
        buttons->csv_load_button.x, buttons->csv_load_button.y,
        buttons->csv_load_button.width, buttons->csv_load_button.height);
    // update hover state for add body button'
    buttons->add_body_button.is_hovered = isMouseInRect(mouse_x, mouse_y,
        buttons->add_body_button.x, buttons->add_body_button.y,
        buttons->add_body_button.width, buttons->add_body_button.height);
    // update hover state for stats button
    buttons->show_stats_button.is_hovered = isMouseInRect(mouse_x, mouse_y,
        buttons->show_stats_button.x, buttons->show_stats_button.y,
        buttons->show_stats_button.width, buttons->show_stats_button.height);
}

// handles mouse button down events (button clicks and drag start)
static void handleMouseButtonDownEvent(SDL_Event* event, window_params_t* wp, body_properties_t** gb, int* num_bodies, spacecraft_properties_t** sc, int* num_craft, button_storage_t* buttons, text_input_dialog_t* dialog, stats_window_t* stats_window) {
    // check if right mouse button or middle mouse button (for dragging)
    if (event->button.button == SDL_BUTTON_RIGHT || event->button.button == SDL_BUTTON_MIDDLE) {
        wp->is_dragging = true;
        wp->drag_start_x = (int)event->button.x;
        wp->drag_start_y = (int)event->button.y;
        wp->drag_origin_x = wp->screen_origin_x;
        wp->drag_origin_y = wp->screen_origin_y;
    }
    // existing button click handling
    else if (buttons->csv_load_button.is_hovered) {
        // reads the CSV file associated with loading orbital bodies
        readCSV("planet_data.csv", gb, num_bodies);
        readSpacecraftCSV("spacecraft_data.csv", sc, num_craft);
        wp->sim_time = 0;
        // reset energy measurement so stats recalculate with new bodies
        stats_window->measured_initial_energy = false;
    }
    else if(buttons->add_body_button.is_hovered) {
        // activate text input dialog
        dialog->active = true;
        dialog->state = INPUT_NAME;
        dialog->input_buffer[0] = '\0';
        SDL_StartTextInput(SDL_GetKeyboardFocus());
    }
    else if(buttons->show_stats_button.is_hovered) {
        // toggle stats display in main window
        stats_window->is_shown = !stats_window->is_shown;
    }
}

// handles mouse button release events
static void handleMouseButtonUpEvent(SDL_Event* event, window_params_t* wp) {
    if (event->button.button == SDL_BUTTON_RIGHT || event->button.button == SDL_BUTTON_MIDDLE) {
        wp->is_dragging = false;
    }
}

// handles mouse wheel events (zooming and speed control)
static void handleMouseWheelEvent(SDL_Event* event, window_params_t* wp, button_storage_t* buttons) {
    int mouse_x = (int)event->wheel.mouse_x;
    int mouse_y = (int)event->wheel.mouse_y;

    // if its in the speed control button area
    if (isMouseInRect(mouse_x, mouse_y,
        buttons->sc_button.x, buttons->sc_button.y,
        buttons->sc_button.width, buttons->sc_button.height)) {

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
static void handleKeyboardEvent(SDL_Event* event, window_params_t* wp, body_properties_t** gb, int* num_bodies, spacecraft_properties_t** sc, int* num_craft) {
    if(event->key.key == SDLK_SPACE) {
        if (wp->sim_running == false) {
            wp->sim_running = true;
        }
        else if (wp->sim_running == true) {
            wp->sim_running = false;
        }
    }
    else if (event->key.key == SDLK_R) {
        resetSim(&wp->sim_time, gb, num_bodies, sc, num_craft);
    }
}

// handles window resize events
static void handleWindowResizeEvent(SDL_Event* event, window_params_t* wp, button_storage_t* buttons) {
    wp->window_size_x = event->window.data1;
    wp->window_size_y = event->window.data2;
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
void runEventCheck(SDL_Event* event, window_params_t* wp, body_properties_t** gb, int* num_bodies, spacecraft_properties_t** sc, int* num_craft, button_storage_t* buttons, text_input_dialog_t* dialog, stats_window_t* stats_window) {
    while (SDL_PollEvent(event)) {
        // if dialog is active, handle text input events
        if (dialog->active) {
            if (handleTextInputDialogEvent(event, wp, dialog, gb, num_bodies)) {
                continue; // skip other event processing when dialog is active
            }
            continue;
        }

        // check if x button is pressed to quit
        if (event->type == SDL_EVENT_QUIT) {
            wp->window_open = false;
        }
        // check if mouse is moving to update hover state
        else if (event->type == SDL_EVENT_MOUSE_MOTION && event->window.windowID == wp->main_window_ID) {
            handleMouseMotionEvent(event, wp, buttons);
        }
        // check if mouse button is clicked
        else if (event->type == SDL_EVENT_MOUSE_BUTTON_DOWN && event->window.windowID == wp->main_window_ID) {
            handleMouseButtonDownEvent(event, wp, gb, num_bodies, sc, num_craft, buttons, dialog, stats_window);
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
            handleKeyboardEvent(event, wp, gb, num_bodies, sc, num_craft);
        }
        // check if window is resized (only for main window)
        else if (event->type == SDL_EVENT_WINDOW_RESIZED &&
                 event->window.windowID == wp->main_window_ID) {
            handleWindowResizeEvent(event, wp, buttons);
        }
    }
}