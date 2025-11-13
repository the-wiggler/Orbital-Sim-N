////////////////////////////////////////////////////////////////////////////////////////////////////
// SIMULATION CALCULATION FUNCTIONS
////////////////////////////////////////////////////////////////////////////////////////////////////
#define SDL_MAIN_HANDLED
#include "calculation_functions.h"
SDL_Color text_color = {255, 255, 255, 255};

// at the end of each sim loop, this function should be run to calculate the changes in
// the force values based on other parameters. for example, using F to find a based on m.
// b is the body that has the force applied to it, whilst b2 is the body applying force to b
void calculateForce(body_properties_t *b, body_properties_t b2) {
    // calculate the distance between the two bodies
    double delta_pos_x = b2.pos_x - b->pos_x;
    double delta_pos_y = b2.pos_y - b->pos_y;
    b->r_from_body = sqrt(delta_pos_x * delta_pos_x + delta_pos_y * delta_pos_y);
    double r = b->r_from_body;

    // calculate the force that b2 applies on b due to gravitation (F = (GMm) / r)
    double total_force = (G * b->mass * b2.mass) / (r * r);
    b->force_x += total_force * (delta_pos_x / r);
    b->force_y += total_force * (delta_pos_y / r);
}

// this calculates the changes of velocity and position based on the force values from before
void updateMotion(body_properties_t *b, double dt) {
    // calculate the acceleration from the force on the object
    b->acc_x = b->force_x / b->mass;
    b->acc_y = b->force_y / b->mass;

    // update the velocity
    b->vel_x = b->vel_x + b->acc_x * dt;
    b->vel_y = b->vel_y + b->acc_y * dt;
    b->vel   = sqrt(b->vel_x * b->vel_x + b->vel_y * b->vel_y);

    // update the position using the new velocity
    b->pos_x = b->pos_x + b->vel_x * dt;
    b->pos_y = b->pos_y + b->vel_y * dt;
}

// transforms spacial coordinates (for example, in meters) to pixel coordinates
void transformCoordinates(body_properties_t *b, window_params_t wp) {
    b->pixel_coordinates_x = wp.screen_origin_x + (int)(b->pos_x / wp.meters_per_pixel);
    b->pixel_coordinates_y = wp.screen_origin_y - (int)(b->pos_y / wp.meters_per_pixel); // this is negative because the SDL origin is in the top left, so positive y is 'down'
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
void drawScaleBar(SDL_Renderer* renderer, double meters_per_pixel, int window_width, int window_height) {
    const int BAR_HEIGHT = 3, MARGIN = 20, BAR_WIDTH_PIXELS = 150;
    
    double distance_meters = BAR_WIDTH_PIXELS * meters_per_pixel;
    
    // determine magnitude and round to nice value
    double magnitude = pow(10.0, floor(log10(distance_meters)));
    double scale_value = round(distance_meters / magnitude) * magnitude;
    
    // format text based on scale
    char scale_text[50];
    sprintf(scale_text, distance_meters >= 1000000 ? "%.0f km" : "%.0f m", 
            distance_meters >= 1000000 ? scale_value / 1000.0 : scale_value);
    
    // draw scale bar
    int bar_x = MARGIN, bar_y = window_height - MARGIN - BAR_HEIGHT;
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


// calculates the size (in pixels) that the planet should appear on the screen based on its mass
int calculateVisualRadius(body_properties_t body, window_params_t wp) {
    int r = (int)(body.radius / wp.meters_per_pixel);
    return r;
}

// thing that calculates changing sim speed
bool isMouseInRect(int mouse_x, int mouse_y, int rect_x, int rect_y, int rect_w, int rect_h) {
    return (mouse_x >= rect_x && mouse_x <= rect_x + rect_w &&
            mouse_y >= rect_y && mouse_y <= rect_y + rect_h);
}

// the speed control box
void drawSpeedControl(SDL_Renderer* renderer, speed_control_t* control, double multiplier) {
    // background color
    SDL_SetRenderDrawColor(renderer, 80, 80, 120, 255);
    
    SDL_FRect bg_rect = {
        (float)control->x, 
        (float)control->y, 
        (float)control->width, 
        (float)control->height
    };
    SDL_RenderFillRect(renderer, &bg_rect);
    
    // text showing current speed
    char speed_text[32];
    snprintf(speed_text, sizeof(speed_text), "Speed: %.2f s/frame", multiplier);
    
    SDL_WriteText(renderer, g_font, speed_text, control->x + 10, control->y + 10, text_color);

}

// the event handling code... checks if events are happening for input and does a task based on that input
void runEventCheck(SDL_Event* event, bool* loop_running_condition, speed_control_t* speed_control, window_params_t* wp, bool* sim_running) {
    while (SDL_PollEvent(event)) {
        // check if x button is pressed to quit
        if (event->type == SDL_EVENT_QUIT) {
            *loop_running_condition = false;
        }
        // check if window is resized
        else if (event->type == SDL_EVENT_WINDOW_RESIZED) {
            // get new window dimensions
            wp->window_size_x = event->window.data1;
            wp->window_size_y = event->window.data2;

            // recalculate screen origin to keep it centered
            wp->screen_origin_x = wp->window_size_x / 2;
            wp->screen_origin_y = wp->window_size_y / 2;

            // recalculate font size based on new window size
            wp->font_size = wp->window_size_x / (wp->window_size_x * 0.05);
        }
        // check if scroll
        else if (event->type == SDL_EVENT_MOUSE_WHEEL) {
            int mouse_x = (int)event->wheel.mouse_x;
            int mouse_y = (int)event->wheel.mouse_y;
            
            // if its in the speed change box
            if (isMouseInRect(mouse_x, mouse_y,
                speed_control->x, speed_control->y,
                speed_control->width, speed_control->height)) {
                
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
                if (*sim_running == false) {
                    *sim_running = true;
                }
                else if (*sim_running == true) {
                    *sim_running = false;
                }
            }
        }

    }
}

// the stats box that shows stats yay
void drawStatsBox(SDL_Renderer* renderer, body_properties_t* bodies, int num_bodies, double sim_time, window_params_t wp) {
    // all calculations for things to go inside the box:
    for (int i = 0; i < num_bodies; i++) {
        char vel_text[32];
        snprintf(vel_text, sizeof(vel_text), "Vel of Body %d: %.1f", i, bodies[i].vel);

        // render text

        SDL_WriteText(renderer, g_font, vel_text, 20, 70+i*20, text_color);
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
    SDL_WriteText(renderer, g_font, time, wp.window_size_x - (wp.window_size_x * 0.2), 15, text_color);
}


// function to add a new body to the system
void addOrbitalBody(body_properties_t** gb, int* num_bodies, double mass, double x_pos, double y_pos, double x_vel, double y_vel) {

    // reallocate memory for the new body
    *gb = (body_properties_t *)realloc(*gb, (*num_bodies + 1) * sizeof(body_properties_t));

    // initialize the new body at index num_bodies
    (*gb)[*num_bodies].mass = mass;
    (*gb)[*num_bodies].pos_x = x_pos;
    (*gb)[*num_bodies].pos_y = y_pos;
    (*gb)[*num_bodies].vel_x = x_vel;
    (*gb)[*num_bodies].vel_y = y_vel;
    
    // calculate the radius based on mass
    (*gb)[*num_bodies].radius = pow(mass, 0.279f);

    // increment the body count
    (*num_bodies)++;
}

// function to remove a body from the system
