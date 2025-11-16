#include "sdl_elements.h"
#include "calculation_functions.h"

SDL_Color text_color = {255, 255, 255, 255};

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
    
    SDL_Surface* text_surface = TTF_RenderText_Solid(font, text, 0, color);
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

// draw a distance scale bar on the sreen
void drawScaleBar(SDL_Renderer* renderer, window_params_t wp) {
    const int BAR_HEIGHT = wp.window_size_y * 0.003;
    const int MARGIN = wp.window_size_x * 0.02;
    const int BAR_WIDTH_PIXELS = wp.window_size_x * 0.15;

    double distance_meters = BAR_WIDTH_PIXELS * wp.meters_per_pixel;
    
    // determine magnitude and round to nice value
    double magnitude = pow(10.0, floor(log10(distance_meters)));
    double scale_value = round(distance_meters / magnitude) * magnitude;
    
    // format text based on scale
    char scale_text[50];
    sprintf(scale_text, distance_meters >= 1000000 ? "%.0f km" : "%.0f m", 
            distance_meters >= 1000000 ? scale_value / 1000.0 : scale_value);
    
    // draw scale bar
    int bar_x = MARGIN; 
    int bar_y = wp.window_size_y - MARGIN - BAR_HEIGHT;
    SDL_SetRenderDrawColor(renderer, 255, 255, 255, 255);
    
    SDL_FRect rects[3] = {
        {bar_x, bar_y, BAR_WIDTH_PIXELS, BAR_HEIGHT},              // main bar
        {bar_x, bar_y - 3, 2, BAR_HEIGHT + 6},                     // left cap
        {bar_x + BAR_WIDTH_PIXELS - 2, bar_y - 3, 2, BAR_HEIGHT + 6}  // right cap
    };
    
    for (int i = 0; i < 3; i++) {
        SDL_RenderFillRect(renderer, &rects[i]);
    }
    
    // Render text
    SDL_WriteText(renderer, g_font, scale_text, bar_x, bar_y, (SDL_Color){255, 255, 255, 255});
}

// thing that calculates changing sim speed
bool isMouseInRect(int mouse_x, int mouse_y, int rect_x, int rect_y, int rect_w, int rect_h) {
    return (mouse_x >= rect_x && mouse_x <= rect_x + rect_w &&
            mouse_y >= rect_y && mouse_y <= rect_y + rect_h);
}

void renderOrbitBodies(SDL_Renderer* renderer, body_properties_t* gb, int num_bodies, window_params_t wp) {
    for (int i = 0; i < num_bodies; i++) {
        // draw bodies
        SDL_RenderFillCircle(renderer, gb[i].pixel_coordinates_x,
                        gb[i].pixel_coordinates_y, 
                        calculateVisualRadius(gb[i], wp));
    }
}

// the stats box that shows stats yay
void drawStatsBox(SDL_Renderer* renderer, body_properties_t* bodies, int num_bodies, double sim_time, window_params_t wp) {
    int margin_x = wp.window_size_x * 0.02;
    int start_y = wp.window_size_y * 0.07;
    int line_height = wp.window_size_y * 0.02;

    // all calculations for things to go inside the box:
    if (bodies != NULL) {
        for (int i = 0; i < num_bodies; i++) {
            char vel_text[32];
            snprintf(vel_text, sizeof(vel_text), "Vel of %s: %.1f", bodies[i].name, bodies[i].vel);

            // render text
            SDL_WriteText(renderer, g_font, vel_text, margin_x, start_y + i * line_height, text_color);
        }
    }

    // show sim time in the top corner
    char time[32];
    if (sim_time < 60) {
        snprintf(time, sizeof(time), "Sim time: %.f s", sim_time); // sim time in seconds
    }
    else if (sim_time > 60 && sim_time < 3600) {
        snprintf(time, sizeof(time), "Sim time: %.f mins", sim_time/60); // sim time in minutes
    }
    else if (sim_time > 3600 && sim_time < 86400) {
        snprintf(time, sizeof(time), "Sim time: %.1f hrs", sim_time/3600); // sim time in hours
    }
    else {
        snprintf(time, sizeof(time), "Sim time: %.1f days", sim_time/86400); // sim time in days
    }
    SDL_WriteText(renderer, g_font, time, wp.window_size_x * 0.75, wp.window_size_y * 0.04, text_color);

    // show paused indication
    if (wp.sim_running) SDL_WriteText(renderer, g_font, "Sim Running...", wp.window_size_x * 0.75, wp.window_size_y * 0.015, text_color);
    else SDL_WriteText(renderer, g_font, "Sim Paused", wp.window_size_x * 0.75, wp.window_size_y * 0.015, text_color);
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
    SDL_WriteText(renderer, g_font, text, button->x + padding_x, button->y + padding_y, text_color);
}

// renders all of the buttons on the screen, this function holds all button drawing logic
void renderUIButtons(SDL_Renderer* renderer, button_storage_t* buttons, window_params_t* wp) {
    // speed control button
    char speed_text[32];
    snprintf(speed_text, sizeof(speed_text), "Speed: %.2f s/frame", wp->time_step);
    renderButton(renderer, &buttons->sc_button, speed_text, *wp);
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// EVENT CHECKING FUNCTION
////////////////////////////////////////////////////////////////////////////////////////////////////
// the event handling code... checks if events are happening for input and does a task based on that input
void runEventCheck(SDL_Event* event, window_params_t* wp, body_properties_t** bodies, int* num_bodies, button_storage_t* buttons) {
    while (SDL_PollEvent(event)) {
        // check if x button is pressed to quit
        if (event->type == SDL_EVENT_QUIT) {
            wp->window_open = false;
        }
        // check if mouse button is clicked
        else if (event->type == SDL_EVENT_MOUSE_BUTTON_DOWN) {

        }
        // check if mouse is moving to update hover state
        else if (event->type == SDL_EVENT_MOUSE_MOTION) {
            int mouse_x = (int)event->motion.x;
            int mouse_y = (int)event->motion.y;

            // update hover state for speed control button
            buttons->sc_button.is_hovered = isMouseInRect(mouse_x, mouse_y,
                buttons->sc_button.x, buttons->sc_button.y,
                buttons->sc_button.width, buttons->sc_button.height);
        }
        // check if scroll
        else if (event->type == SDL_EVENT_MOUSE_WHEEL) {
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
        // check if keyboard key is pressed
        else if (event->type == SDL_EVENT_KEY_DOWN) {
            if(event->key.key == SDLK_SPACE) {
                if (wp->sim_running == false) {
                    wp->sim_running = true;
                }
                else if (wp->sim_running == true) {
                    wp->sim_running = false;
                }
            }
            else if (event->key.key == SDLK_R) {
                resetSim(&wp->sim_time, bodies, num_bodies);
            }
        }
        // check if window is resized
        else if (event->type == SDL_EVENT_WINDOW_RESIZED) {
            wp->window_size_x = event->window.data1;
            wp->window_size_y = event->window.data2;
            wp->screen_origin_x = wp->window_size_x / 2;
            wp->screen_origin_y = wp->window_size_y / 2;

            // update font size proportionally
            wp->font_size = wp->window_size_x / 75;
            if (g_font) TTF_CloseFont(g_font);
            g_font = TTF_OpenFont("CascadiaCode.ttf", wp->font_size);

            // update speed control button position and size
            buttons->sc_button.width = wp->window_size_x * 0.25;
            buttons->sc_button.height = wp->window_size_y * 0.04;
            buttons->sc_button.x = wp->window_size_x * 0.01;
            buttons->sc_button.y = wp->window_size_y * 0.007;
        }

    }
}