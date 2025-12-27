#include "../gui/SDL_engine.h"
#include "../globals.h"
#include "../sim/bodies.h"
#include <SDL3/SDL.h>
#include <GL/glew.h>
#include <GL/gl.h>

// display error message using SDL dialog
void displayError(const char* title, const char* message) {
    SDL_ShowSimpleMessageBox(SDL_MESSAGEBOX_ERROR, title, message, NULL);
}
// initialize the window parameters
window_params_t init_window_params() {
    window_params_t wp = {0};

    // gets screen width and height parameters from computer
    const SDL_DisplayMode *mode = SDL_GetCurrentDisplayMode(SDL_GetPrimaryDisplay());

    wp.time_step = 1;
    // sets the default window size scaled based on the user's screen size
    wp.window_size_x = (float)mode->w * (2.0f/3.0f);
    wp.window_size_y = (float)mode->h * (2.0f/3.0f);

    // initialize 3D camera
    wp.camera_pos.x = 2.0f; wp.camera_pos.y = 2.0f; wp.camera_pos.z = 3.0f;
    wp.zoom = 1.5f;

    wp.meters_per_pixel = 100000.0;

    wp.font_size = (float)wp.window_size_x / 50;
    wp.window_open = true;
    wp.sim_running = true;
    wp.data_logging_enabled = false;
    wp.sim_time = 0;

    wp.is_dragging = false;

    wp.main_view_shown = true;
    wp.craft_view_shown = false;

    return wp;
}

SDL_GL_init_t init_SDL_OPENGL_window(const char* title, int width, int height, Uint32* outWindowID) {
    SDL_GL_init_t result = {0};

    // set OpenGL attributes
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MAJOR_VERSION, 3);
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MINOR_VERSION, 3);
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_PROFILE_MASK, SDL_GL_CONTEXT_PROFILE_CORE);
    SDL_GL_SetAttribute(SDL_GL_DOUBLEBUFFER, 1);
    SDL_GL_SetAttribute(SDL_GL_DEPTH_SIZE, 24);
    SDL_GL_SetAttribute(SDL_GL_STENCIL_SIZE, 8);

    // create SDL window
    result.window = SDL_CreateWindow(title, width, height,
        SDL_WINDOW_OPENGL | SDL_WINDOW_RESIZABLE);

    // store window ID
    if (outWindowID) {
        *outWindowID = SDL_GetWindowID(result.window);
    }

    // initialize OpenGL context and GLEW
    result.glContext = SDL_GL_CreateContext(result.window);
    glewExperimental = GL_TRUE;
    glewInit();

    // enable VSync
    SDL_GL_SetSwapInterval(1);

    printf("OpenGL version: %s\n", glGetString(GL_VERSION)); // yes this gives a warning, its supposed to be like this
    printf("GLEW version: %s\n", glewGetString(GLEW_VERSION));

    return result;
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

    // viewport dragging
    if (wp->is_dragging) {
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

    if (event->wheel.y > 0) {
        wp->is_zooming = true;
        wp->zoom *= 1.1f;
    } else if (event->wheel.y < 0) {
        wp->is_zooming = true;
        wp->zoom /= 1.1f;
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
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// MAIN EVENT CHECKING FUNCTION
////////////////////////////////////////////////////////////////////////////////////////////////////
// the event handling code... checks if events are happening for input and does a task based on that input
void runEventCheck(SDL_Event* event, sim_properties_t* sim) {
    window_params_t* wp = &sim->wp;

    wp->is_zooming = false; // reset zooming checks
    wp->is_zooming_in = false;
    wp->is_zooming_out = false;

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
