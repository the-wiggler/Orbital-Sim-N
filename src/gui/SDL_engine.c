#include "../gui/SDL_engine.h"
#include "../globals.h"
#include "../sim/bodies.h"

// display error message using SDL dialog
void displayError(const char* title, const char* message) {
    SDL_ShowSimpleMessageBox(SDL_MESSAGEBOX_ERROR, title, message, NULL);
}

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

// thing that calculates changing sim speed
bool isMouseInRect(const int mouse_x, const int mouse_y, const int rect_x, const int rect_y, const int rect_w, const int rect_h) {
    return (mouse_x >= rect_x && mouse_x <= rect_x + rect_w &&
            mouse_y >= rect_y && mouse_y <= rect_y + rect_h);
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// EVENT HANDLING HELPER FUNCTIONS
////////////////////////////////////////////////////////////////////////////////////////////////////
// handles mouse motion events (dragging and button hover states)
static void handleMouseMotionEvent(const SDL_Event* event, sim_properties_t* sim) {
    window_params_t* wp = &sim->wp;

    int mouse_x = (int)event->motion.x;
    int mouse_y = (int)event->motion.y;

    // viewport dragging
    if (wp->is_dragging) {
        int delta_x = mouse_x - (int)wp->drag_start_x;
        int delta_y = mouse_y - (int)wp->drag_start_y;
        wp->screen_origin_x = wp->drag_origin_x + (float)delta_x;
        wp->screen_origin_y = wp->drag_origin_y + (float)delta_y;
    }
}

// handles mouse button down events (button clicks and drag start)
static void handleMouseButtonDownEvent(const SDL_Event* event, sim_properties_t* sim) {
    window_params_t* wp = &sim->wp;
    body_properties_t* gb = &sim->gb;
    spacecraft_properties_t* sc = &sim->gs;

    // check if right mouse button or middle mouse button (for dragging)
    if (event->button.button == SDL_BUTTON_RIGHT || event->button.button == SDL_BUTTON_MIDDLE) {
        wp->is_dragging = true;
        wp->drag_start_x = event->button.x;
        wp->drag_start_y = event->button.y;
        wp->drag_origin_x = wp->screen_origin_x;
        wp->drag_origin_y = wp->screen_origin_y;
    }
}

// handles mouse button release events
static void handleMouseButtonUpEvent(const SDL_Event* event, sim_properties_t* sim) {
    window_params_t* wp = &sim->wp;

    if (event->button.button == SDL_BUTTON_RIGHT || event->button.button == SDL_BUTTON_MIDDLE) {
        wp->is_dragging = false;
    }
}

// handles mouse wheel events (zooming and speed control)
static void handleMouseWheelEvent(const SDL_Event* event, sim_properties_t* sim) {
    window_params_t* wp = &sim->wp;

    int mouse_x = (int)event->wheel.mouse_x;
    int mouse_y = (int)event->wheel.mouse_y;

    if (event->wheel.y > 0) {
        wp->is_zooming = true;
        wp->meters_per_pixel *= 1.05;
    } else if (event->wheel.y < 0) {
        wp->is_zooming = true;
        wp->meters_per_pixel /= 1.05;
    }
}

// handles keyboard events (pause/play, reset)
static void handleKeyboardEvent(const SDL_Event* event, sim_properties_t* sim) {
    window_params_t* wp = &sim->wp;

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
    else if (event->key.key == SDLK_D) {
        if (wp->data_logging_enabled) {
            wp->data_logging_enabled = false;
        }
        else if (!wp->data_logging_enabled) {
            wp->data_logging_enabled = true;
        }
    }
}

// handles window resize events
static void handleWindowResizeEvent(const SDL_Event* event, sim_properties_t* sim) {
    window_params_t* wp = &sim->wp;

    wp->window_size_x = (float)event->window.data1;
    wp->window_size_y = (float)event->window.data2;
    wp->screen_origin_x = wp->window_size_x / 2;
    wp->screen_origin_y = wp->window_size_y / 2;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// MAIN EVENT CHECKING FUNCTION
////////////////////////////////////////////////////////////////////////////////////////////////////
// the event handling code... checks if events are happening for input and does a task based on that input
void runEventCheck(SDL_Event* event, sim_properties_t* sim) {
    window_params_t* wp = &sim->wp;

    while (SDL_PollEvent(event)) {
        // check if x button is pressed to quit
        if (event->type == SDL_EVENT_QUIT) {
            wp->reset_sim = true;
            wp->window_open = false;
            wp->sim_running = false;
        }
        // check if mouse is moving to update hover state
        else if (event->type == SDL_EVENT_MOUSE_MOTION && event->window.windowID == wp->main_window_ID) {
            handleMouseMotionEvent(event, sim);
        }
        // check if mouse button is clicked
        else if (event->type == SDL_EVENT_MOUSE_BUTTON_DOWN && event->window.windowID == wp->main_window_ID) {
            handleMouseButtonDownEvent(event, sim);
        }
        // check if mouse button is released
        else if (event->type == SDL_EVENT_MOUSE_BUTTON_UP && event->window.windowID == wp->main_window_ID) {
            handleMouseButtonUpEvent(event, sim);
        }
        // check if scroll
        else if (event->type == SDL_EVENT_MOUSE_WHEEL && event->window.windowID == wp->main_window_ID) {
            handleMouseWheelEvent(event, sim);
        }
        // check if keyboard key is pressed
        else if (event->type == SDL_EVENT_KEY_DOWN && event->window.windowID == wp->main_window_ID) {
            handleKeyboardEvent(event, sim);
        }
        // check if window is resized (only for main window)
        else if (event->type == SDL_EVENT_WINDOW_RESIZED &&
                 event->window.windowID == wp->main_window_ID) {
            handleWindowResizeEvent(event, sim);
        }
    }
}
