#include "../gui/SDL_engine.h"

#include <stdlib.h>

#include "../gui/GL_renderer.h"
#include <SDL3/SDL.h>
#include <GL/glew.h>

#include "../globals.h"
#include "../utility/json_loader.h"
#ifdef __APPLE__
#include <OpenGL/gl.h>
#else
#include <GL/gl.h>
#endif
#include "../math/matrix.h"
#include <string.h>

// display error message using SDL dialog
void displayError(const char* title, const char* message) {
    SDL_ShowSimpleMessageBox(SDL_MESSAGEBOX_ERROR, title, message, NULL);
}
// initialize the window parameters
window_params_t init_window_params(void) {
    window_params_t wp = {0};

    // gets screen width and height parameters from computer
    const SDL_DisplayMode *mode = SDL_GetCurrentDisplayMode(SDL_GetPrimaryDisplay());

    wp.time_step = 1;
    // sets the default window size scaled based on the user's screen size
    wp.window_size_x = (float)mode->w * (2.0f/3.0f);
    wp.window_size_y = (float)mode->h * (2.0f/3.0f);

    // initialize 3D camera
    wp.cam_yaw = PI_OVER_4_f;
    wp.cam_pitch = PI_OVER_6_f;

    // compute initial camera position from angles
    const float cos_pitch = cosf(wp.cam_pitch);
    wp.camera_pos.x = cos_pitch * cosf(wp.cam_yaw);
    wp.camera_pos.y = cos_pitch * sinf(wp.cam_yaw);
    wp.camera_pos.z = sinf(wp.cam_pitch);

    wp.zoom = 1.5f;

    wp.meters_per_pixel = 100000.0;

    wp.window_open = true;
    wp.sim_running = false;  // start paused until setup is complete
    wp.data_logging_enabled = false;
    wp.sim_time = 0;

    wp.is_dragging_orbit_view = false;
    wp.is_dragging_translation_view = false;

    // initial visual defaults
    wp.draw_lines_between_bodies = false;
    wp.draw_inclination_height = true;
    wp.draw_planet_path = true;
    wp.draw_craft_path = true;
    wp.draw_planet_SOI = true;

    wp.cam_target.x = 0; wp.cam_target.y = 0; wp.cam_target.z = 0;

    wp.frame_counter = 0;

    return wp;
}

console_t init_console(const window_params_t wp) {
    console_t console = {0};

    // init text input
    console.cmd_text_box[0] = '\0';
    console.cmd_text_box_length = 0;
    console.cmd_pos_x = 0.02f * wp.window_size_x;
    console.cmd_pos_y = wp.window_size_y - 0.1f * wp.window_size_y;

    // init console log box
    console.log[0] = '\0';
    console.log_pos_x = console.cmd_pos_x;
    console.log_pos_y = console.cmd_pos_y + 0.05f * wp.window_size_y;

    return console;
}

SDL_GL_init_t init_SDL_OPENGL_window(const char* title, const int width, const int height, Uint32* outWindowID) {
    SDL_GL_init_t result = {0};

    // set OpenGL attributes
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MAJOR_VERSION, 3);
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MINOR_VERSION, 3);
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_PROFILE_MASK, SDL_GL_CONTEXT_PROFILE_CORE);
#ifdef __APPLE__
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_FLAGS, SDL_GL_CONTEXT_FORWARD_COMPATIBLE_FLAG);
#endif
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
    if (!result.glContext) {
        fprintf(stderr, "Failed to create OpenGL context: %s\n", SDL_GetError());
        return result;
    }

    glewExperimental = GL_TRUE;
    const GLenum glewError = glewInit();
    if (glewError != GLEW_OK) {
        fprintf(stderr, "Error initializing GLEW: %s\n", (const char*)glewGetErrorString(glewError));
        return result;
    }

    // enable VSync
    SDL_GL_SetSwapInterval(1);

    printf("OpenGL version: %s\n", (const char*)glGetString(GL_VERSION));
    printf("GLEW version: %s\n", (const char*)glewGetString(GLEW_VERSION));

    // enable text input
    SDL_StartTextInput(result.window);

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

    // orbit camera
    if (wp->is_dragging_orbit_view) {
        // calculate mouse movement delta
        const float delta_x = event->motion.x - wp->drag_last.x;
        const float delta_y = event->motion.y - wp->drag_last.y;

        // convert mouse movement to rotation angles
        const float sensitivity = 0.005f;
        wp->cam_yaw -= delta_x * sensitivity;
        wp->cam_pitch += delta_y * sensitivity;

        const float pitch_limit = PI_OVER_2_f - 0.01f;
        if (wp->cam_pitch > pitch_limit) wp->cam_pitch = pitch_limit;
        if (wp->cam_pitch < -pitch_limit) wp->cam_pitch = -pitch_limit;

        const float cos_pitch = cosf(wp->cam_pitch);
        wp->camera_pos.x = cos_pitch * cosf(wp->cam_yaw);
        wp->camera_pos.y = cos_pitch * sinf(wp->cam_yaw);
        wp->camera_pos.z = sinf(wp->cam_pitch);

        // update last mouse position for next frame
        wp->drag_last.x = event->motion.x;
        wp->drag_last.y = event->motion.y;
    }

    // translation view
    if (wp->is_dragging_translation_view) {
        const float delta_x = event->motion.x - wp->drag_last.x;
        const float delta_y = event->motion.y - wp->drag_last.y;

        // compute camera right vector projected onto XY plane
        const vec3_f forward = {-wp->camera_pos.x, -wp->camera_pos.y, 0.0f};
        const vec3_f world_up = {0.0f, 0.0f, 1.0f};
        const vec3_f right = vec3_f_normalize(vec3_f_cross(forward, world_up));
        // "up" in screen space maps to forward direction in XY plane
        const vec3_f forward_xy = vec3_f_normalize((vec3_f){-wp->camera_pos.x, -wp->camera_pos.y, 0.0f});

        // scale movement by zoom level
        const float sensitivity = 0.001f * wp->zoom;

        // move target in XY plane
        wp->cam_target.x -= right.x * delta_x * sensitivity - forward_xy.x * delta_y * sensitivity;
        wp->cam_target.y -= right.y * delta_x * sensitivity - forward_xy.y * delta_y * sensitivity;

        wp->drag_last.x = event->motion.x;
        wp->drag_last.y = event->motion.y;
    }
}

static void handleMouseButtonDownEvent(const SDL_Event* event, sim_properties_t* sim) {
    window_params_t* wp = &sim->wp;
    // body_properties_t* gb = &sim->gb;
    // spacecraft_properties_t* sc = &sim->gs;

    // check if right mouse button (for camera orbit view)
    if (event->button.button == SDL_BUTTON_RIGHT) {
        wp->is_dragging_orbit_view = true;
        wp->drag_last.x = event->button.x;
        wp->drag_last.y = event->button.y;
    }
    else if (event->button.button == SDL_BUTTON_MIDDLE) {
        wp->is_dragging_translation_view = true;
        wp->drag_last.x = event->button.x;
        wp->drag_last.y = event->button.y;
    }
}

// handles mouse button release events
static void handleMouseButtonUpEvent(const SDL_Event* event, sim_properties_t* sim) {
    window_params_t* wp = &sim->wp;

    if (event->button.button == SDL_BUTTON_RIGHT) wp->is_dragging_orbit_view = false;
    else if (event->button.button == SDL_BUTTON_MIDDLE) wp->is_dragging_translation_view = false;
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

static void parseRunCommands(char* cmd, sim_properties_t* sim) {
    console_t* console = &sim->console;

    if (strncmp(cmd, "step ", 4) == 0) {
        char* argument = cmd + 4;
        sim->wp.time_step = strtod(argument, &argument);

        sprintf(console->log, "step set to %f", sim->wp.time_step);
    }
    else if (strcmp(cmd, "pause") == 0 || strcmp(cmd, "p") == 0) {
        sim->wp.sim_running = false;
        sprintf(console->log, "sim paused");
    }
    else if (strcmp(cmd, "resume") == 0 || strcmp(cmd, "r") == 0) {
        sim->wp.sim_running = true;
        sprintf(console->log, "sim resumed");
    }
    else if (strcmp(cmd, "load") == 0) {
        if (sim->gb.count == 0) {
            sim->wp.sim_running = false; // pauses before loading
            readSimulationJSON("simulation_data.json", &sim->gb, &sim->gs);
            sprintf(console->log, "%d planets and %d craft loaded from json file", sim->gb.count, sim->gs.count);
        }
        else sprintf(console->log, "Warning: system already loaded, reset before loading another");
    }
    else if (strcmp(cmd, "reset") == 0) {
        sim->wp.reset_sim = true;
        sprintf(console->log, "sim reset");
    }
    else if (strncmp(cmd, "enable ", 7) == 0) {
        char* argument = cmd + 7;
        if (strcmp(argument, "guidance-lines") == 0) {
            sim->wp.draw_lines_between_bodies = true;
            sprintf(console->log, "enabled guidance lines");
        }
        else if (strcmp(argument, "export") == 0) {
            sim->wp.data_logging_enabled = true;
            sprintf(console->log, "enabled telemetry export");
        }
        else sprintf(console->log, "unknown argument after command: %s", argument);
    }
    else if (strncmp(cmd, "disable ", 8) == 0) {
        char* argument = cmd + 8;
        if (strcmp(argument, "guidance-lines") == 0) {
            sim->wp.draw_lines_between_bodies = false;
            sprintf(console->log, "disabled guidance lines");
        }
        else if (strcmp(argument, "export") == 0) {
            sim->wp.data_logging_enabled = false;
            sprintf(console->log, "disabled telemetry export");
        }
        else sprintf(console->log, "unknown argument after disable: %s", argument);
    }
    else {
        sprintf(console->log, "unknown command: %s", cmd);
    }

}

// handles keyboard events
static void handleKeyboardEvent(const SDL_Event* event, sim_properties_t* sim) {
    console_t* console = &sim->console;
    //window_params_t* wp = &sim->wp;

    if (event->key.key == SDLK_BACKSPACE && console->cmd_text_box_length > 0) {
        console->cmd_text_box_length -= 1;
        console->cmd_text_box[console->cmd_text_box_length] = '\0';
    }
    else if (event->key.key == SDLK_RETURN || event->key.key == SDLK_KP_ENTER) {
        // clear log because a new command will show a new log message!
        console->log[0] = '\0';
        // if the enter key is pressed, then the command should be queued!
        parseRunCommands(console->cmd_text_box, sim);
        console->cmd_text_box[0] = '\0';
        console->cmd_text_box_length = 0;
    }
}

// handles text input events
static void handleTextInputEvent(const SDL_Event* event, sim_properties_t* sim) {
    console_t* console = &sim->console;

    const size_t text_len = strlen(event->text.text);
    if (console->cmd_text_box_length + text_len < 255) {
        strcat(console->cmd_text_box, event->text.text);
        console->cmd_text_box_length += (int)text_len;
    }
}

// handles window resize events
static void handleWindowResizeEvent(const SDL_Event* event, sim_properties_t* sim) {
    window_params_t* wp = &sim->wp;
    console_t* console = &sim->console;

    wp->window_size_x = (float)event->window.data1;
    wp->window_size_y = (float)event->window.data2;

    console->cmd_pos_x = 0.02f * wp->window_size_x;
    console->cmd_pos_y = wp->window_size_y - 0.1f * wp->window_size_y;

    console->log_pos_x = console->cmd_pos_x;
    console->log_pos_y = console->cmd_pos_y + 0.05f * wp->window_size_y;
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
        // check if text input
        else if (event->type == SDL_EVENT_TEXT_INPUT && event->window.windowID == wp->main_window_ID) {
            handleTextInputEvent(event, sim);
        }
        // check if window is resized
        else if (event->type == SDL_EVENT_WINDOW_RESIZED &&
                 event->window.windowID == wp->main_window_ID) {
            handleWindowResizeEvent(event, sim);
        }
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// CONSOLE RENDERING
////////////////////////////////////////////////////////////////////////////////////////////////////
void renderCMDWindow(sim_properties_t* sim, font_t* font) {
    console_t* console = &sim->console;
    const window_params_t* wp = &sim->wp;

    // display what's in the text box
    addText(font, console->cmd_pos_x, console->cmd_pos_y, console->cmd_text_box, 1.0f);

    // blinking cursor
    if ((wp->frame_counter / 30) % 2 == 0) {
        char cursor_text[256];
        snprintf(cursor_text, sizeof(cursor_text), "%s_", console->cmd_text_box);
        addText(font, console->cmd_pos_x, console->cmd_pos_y, cursor_text, 1.0f);
    }

    // display whatever text is in the log box (resets when enter is pressed)
    addText(font, console->log_pos_x, console->log_pos_y, console->log, 0.8f);
}
