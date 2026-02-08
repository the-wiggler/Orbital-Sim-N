#include <GL/glew.h>
#ifdef __APPLE__
#include <OpenGL/gl.h>
#else
#include <GL/gl.h>
#endif

#include "../gui/SDL_engine.h"

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <SDL3/SDL_messagebox.h>
#include <SDL3/SDL_video.h>
#include <SDL3/SDL_events.h>
#include <SDL3/SDL_keyboard.h>
#include <SDL3/SDL_error.h>
#include <SDL3/SDL_mouse.h>
#include <SDL3/SDL_keycode.h>
#include <SDL3/SDL_hints.h>
#include <SDL3/SDL_init.h>

#include "../globals.h"
#include "../utility/json_loader.h"
#include "../sim/bodies.h"
#include "../gui/GL_renderer.h"
#include "../math/matrix.h"
#include "../types.h"
#include "../sim/spacecraft.h"

// display error message using SDL dialog
void displayError(const char* title, const char* message) {
    SDL_ShowSimpleMessageBox(SDL_MESSAGEBOX_ERROR, title, message, NULL);
}
// initialize the window parameters
window_params_t init_window_params(void) {
    window_params_t window_params = {0};

    // gets screen width and height parameters from computer
    const SDL_DisplayMode *mode = SDL_GetCurrentDisplayMode(SDL_GetPrimaryDisplay());

    window_params.time_step = 1;
    // sets the default window size scaled based on the user's screen size
    window_params.window_size_x = (float)mode->w * (2.0F/3.0F);
    window_params.window_size_y = (float)mode->h * (2.0F/3.0F);

    // initialize 3D camera
    window_params.cam_yaw = PI_OVER_4_f;
    window_params.cam_pitch = PI_OVER_6_f;

    // compute initial camera position from angles
    const float cos_pitch = cosf(window_params.cam_pitch);
    window_params.camera_pos.x = cos_pitch * cosf(window_params.cam_yaw);
    window_params.camera_pos.y = cos_pitch * sinf(window_params.cam_yaw);
    window_params.camera_pos.z = sinf(window_params.cam_pitch);

    window_params.zoom = 1.5F;

    window_params.meters_per_pixel = 100000.0;

    window_params.window_open = true;
    window_params.sim_running = false;  // start paused until setup is complete
    window_params.data_logging_enabled = false;
    window_params.sim_time = 0;

    window_params.is_dragging_orbit_view = false;
    window_params.is_dragging_translation_view = false;

    // initial visual defaults
    window_params.draw_lines_between_bodies = false;
    window_params.draw_inclination_height = true;
    window_params.draw_planet_path = true;
    window_params.draw_craft_path = true;
    window_params.draw_planet_SOI = true;

    window_params.cam_target.x = 0; window_params.cam_target.y = 0; window_params.cam_target.z = 0;

    window_params.frame_counter = 0;

    return window_params;
}

console_t init_console(const window_params_t window_params) {
    console_t console = {0};

    // init text input
    console.cmd_text_box[0] = '\0';
    console.cmd_text_box_length = 0;
    console.cmd_pos_x = 0.02F * window_params.window_size_x;
    console.cmd_pos_y = window_params.window_size_y - (0.1F * window_params.window_size_y);

    // init console log box
    console.log[0] = '\0';
    console.log_pos_x = console.cmd_pos_x;
    console.log_pos_y = console.cmd_pos_y + (0.05F * window_params.window_size_y);

    return console;
}

SDL_GL_init_t init_SDL_OPENGL_window(const char* title, sim_properties_t* sim) {
    SDL_GL_init_t result = {0};

#ifdef __linux__
    // force X11 on Linux (fixes SDL text input issues on wayland)
    SDL_SetHint(SDL_HINT_VIDEO_DRIVER, "x11");
#endif

    // initialize SDL
    SDL_Init(SDL_INIT_VIDEO);

    // window parameters & command prompt init
    sim->window_params = init_window_params();
    sim->console = init_console(sim->window_params);

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
    result.window = SDL_CreateWindow(title, (int)sim->window_params.window_size_x, (int)sim->window_params.window_size_y, SDL_WINDOW_OPENGL | SDL_WINDOW_RESIZABLE);

    // store window ID
    sim->window_params.main_window_ID = SDL_GetWindowID(result.window);

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
    window_params_t* window_params = &sim->window_params;

    // orbit camera
    if (window_params->is_dragging_orbit_view) {
        // calculate mouse movement delta
        const float delta_x = event->motion.x - window_params->drag_last.x;
        const float delta_y = event->motion.y - window_params->drag_last.y;

        // convert mouse movement to rotation angles
        const float sensitivity = 0.005F;
        window_params->cam_yaw -= delta_x * sensitivity;
        window_params->cam_pitch += delta_y * sensitivity;

        const float pitch_limit = PI_OVER_2_f - 0.01F;
        if (window_params->cam_pitch > pitch_limit) {
            window_params->cam_pitch = pitch_limit;
        }
        if (window_params->cam_pitch < -pitch_limit) {
            window_params->cam_pitch = -pitch_limit;
        }

        const float cos_pitch = cosf(window_params->cam_pitch);
        window_params->camera_pos.x = cos_pitch * cosf(window_params->cam_yaw);
        window_params->camera_pos.y = cos_pitch * sinf(window_params->cam_yaw);
        window_params->camera_pos.z = sinf(window_params->cam_pitch);

        // update last mouse position for next frame
        window_params->drag_last.x = event->motion.x;
        window_params->drag_last.y = event->motion.y;
    }

    // translation view
    if (window_params->is_dragging_translation_view) {
        const float delta_x = event->motion.x - window_params->drag_last.x;
        const float delta_y = event->motion.y - window_params->drag_last.y;

        // compute camera right vector projected onto XY plane
        const vec3_f forward = {-window_params->camera_pos.x, -window_params->camera_pos.y, 0.0F};
        const vec3_f world_up = {0.0F, 0.0F, 1.0F};
        const vec3_f right = vec3_f_normalize(vec3_f_cross(forward, world_up));
        // "up" in screen space maps to forward direction in XY plane
        const vec3_f forward_xy = vec3_f_normalize((vec3_f){-window_params->camera_pos.x, -window_params->camera_pos.y, 0.0F});

        // scale movement by zoom level
        const float sensitivity = 0.001F * window_params->zoom;

        // move target in XY plane
        window_params->cam_target.x -= (right.x * delta_x * sensitivity) - (forward_xy.x * delta_y * sensitivity);
        window_params->cam_target.y -= (right.y * delta_x * sensitivity) - (forward_xy.y * delta_y * sensitivity);

        window_params->drag_last.x = event->motion.x;
        window_params->drag_last.y = event->motion.y;
    }
}

static void handleMouseButtonDownEvent(const SDL_Event* event, sim_properties_t* sim) {
    window_params_t* window_params = &sim->window_params;
    // body_properties_t* gb = &sim->global_bodies;
    // spacecraft_properties_t* sc = &sim->global_spacecraft;

    // check if right mouse button (for camera orbit view)
    if (event->button.button == SDL_BUTTON_RIGHT) {
        window_params->is_dragging_orbit_view = true;
        window_params->drag_last.x = event->button.x;
        window_params->drag_last.y = event->button.y;
    }
    else if (event->button.button == SDL_BUTTON_MIDDLE) {
        window_params->is_dragging_translation_view = true;
        window_params->drag_last.x = event->button.x;
        window_params->drag_last.y = event->button.y;
    }
}

// handles mouse button release events
static void handleMouseButtonUpEvent(const SDL_Event* event, sim_properties_t* sim) {
    window_params_t* window_params = &sim->window_params;

    if (event->button.button == SDL_BUTTON_RIGHT) {
        window_params->is_dragging_orbit_view = false;
    } else if (event->button.button == SDL_BUTTON_MIDDLE) {
        window_params->is_dragging_translation_view = false;
    }
}

// handles mouse wheel events (zooming and speed control)
static void handleMouseWheelEvent(const SDL_Event* event, sim_properties_t* sim) {
    window_params_t* window_params = &sim->window_params;

    if (event->wheel.y > 0) {
        window_params->is_zooming = true;
        window_params->zoom *= 1.1F;
    } else if (event->wheel.y < 0) {
        window_params->is_zooming = true;
        window_params->zoom /= 1.1F;
    }
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
static void parseRunCommands(char* cmd, sim_properties_t* sim, binary_filenames_t* filenames) {
    console_t* console = &sim->console;

    if (strncmp(cmd, "step ", 4) == 0) {
        char* argument = cmd + 4;
        sim->window_params.time_step = strtod(argument, &argument);
        // recompute rotations when timestep changes
        body_precomputeRotations(&sim->global_bodies, sim->window_params.time_step);

        snprintf(console->log, sizeof(console->log), "step set to %f", sim->window_params.time_step);
    }
    else if (strcmp(cmd, "pause") == 0 || strcmp(cmd, "p") == 0) {
        sim->window_params.sim_running = false;
        snprintf(console->log, sizeof(console->log), "sim paused");
    }
    else if (strcmp(cmd, "resume") == 0 || strcmp(cmd, "r") == 0) {
        sim->window_params.sim_running = true;
        snprintf(console->log, sizeof(console->log), "sim resumed");
    }
    else if (strcmp(cmd, "load") == 0) {
        if (sim->global_bodies.count == 0) {
            sim->window_params.sim_running = false; // pauses before loading
            readSimulationJSON("simulation_data.json", &sim->global_bodies, &sim->global_spacecraft);
            // precompute rotations for all bodies based on current timestep
            body_precomputeRotations(&sim->global_bodies, sim->window_params.time_step);
            snprintf(console->log, sizeof(console->log), "%d planets and %d craft loaded from json file", sim->global_bodies.count, sim->global_spacecraft.count);
        } else {
            snprintf(console->log, sizeof(console->log), "Warning: system already loaded, reset before loading another");
        }
    }
    else if (strcmp(cmd, "reset") == 0) {
        sim->window_params.reset_sim = true;
        snprintf(console->log, sizeof(console->log), "sim reset");
    }
    else if (strncmp(cmd, "enable ", 7) == 0) {
        char* argument = cmd + 7;
        if (strcmp(argument, "guidance-lines") == 0) {
            sim->window_params.draw_lines_between_bodies = true;
            snprintf(console->log, sizeof(console->log), "enabled guidance lines");
        }
        else if (strcmp(argument, "export") == 0) {
            sim->window_params.data_logging_enabled = true;
            snprintf(console->log, sizeof(console->log), "enabled telemetry export");
        } else {
            snprintf(console->log, sizeof(console->log), "unknown argument after command: %s", argument);
        }
    }
    else if (strncmp(cmd, "disable ", 8) == 0) {
        char* argument = cmd + 8;
        if (strcmp(argument, "guidance-lines") == 0) {
            sim->window_params.draw_lines_between_bodies = false;
            snprintf(console->log, sizeof(console->log), "disabled guidance lines");
        }
        else if (strcmp(argument, "export") == 0) {
            sim->window_params.data_logging_enabled = false;
            snprintf(console->log, sizeof(console->log), "disabled telemetry export");
        }
        else {
            snprintf(console->log, sizeof(console->log), "unknown argument after disable: %s", argument);
        }
    }
    else if (strncmp(cmd, "sample-period ", 14) == 0) {
        char* argument = cmd + 14;
        filenames->csv_update_period = strtod(argument, &argument);
        snprintf(console->log, sizeof(console->log), "CSV sample period changed to %fs", filenames->csv_update_period);
    }
    else if (strncmp(cmd, "auto ", 5) == 0) {
        char* argument = cmd + 5;
        const int craft_idx = findSpacecraftID(&sim->global_spacecraft, argument);
        if (craft_idx != -1) {
            spacecraft_t* craft = &sim->global_spacecraft.spacecraft[craft_idx];
            craft->burn_properties[craft->num_burns] = craft_autoDeltaVOptimization(sim, craft_idx);
            craft->num_burns++;
            snprintf(console->log, sizeof(console->log), "Optimal target burn created for %s", argument);
        } else {
            snprintf(console->log, sizeof(console->log), "unknown spacecraft: %s", argument);
        }
    }
    else {
        snprintf(console->log, sizeof(console->log), "unknown command: %s", cmd);
    }

}

// handles keyboard events
static void handleKeyboardEvent(const SDL_Event* event, sim_properties_t* sim, binary_filenames_t* filenames) {
    console_t* console = &sim->console;
    //window_params_t* window_params = &sim->wp;

    if (event->key.key == SDLK_BACKSPACE && console->cmd_text_box_length > 0) {
        console->cmd_text_box_length -= 1;
        console->cmd_text_box[console->cmd_text_box_length] = '\0';
    }
    else if (event->key.key == SDLK_RETURN || event->key.key == SDLK_KP_ENTER) {
        // clear log because a new command will show a new log message!
        console->log[0] = '\0';
        // if the enter key is pressed, then the command should be queued!
        parseRunCommands(console->cmd_text_box, sim, filenames);
        console->cmd_text_box[0] = '\0';
        console->cmd_text_box_length = 0;
    }
}

// handles text input events
static void handleTextInputEvent(const SDL_Event* event, sim_properties_t* sim) {
    console_t* console = &sim->console;

    const size_t text_len = strlen(event->text.text);
    if (console->cmd_text_box_length + text_len < 255) {
        #ifdef _WIN32
        strcat_s(console->cmd_text_box, sizeof(console->cmd_text_box), event->text.text);
        #else
        strncat(console->cmd_text_box, event->text.text, sizeof(console->cmd_text_box) - console->cmd_text_box_length - 1);
        #endif
        console->cmd_text_box_length += (int)text_len;
    }
}

// handles window resize events
static void handleWindowResizeEvent(const SDL_Event* event, sim_properties_t* sim) {
    window_params_t* window_params = &sim->window_params;
    console_t* console = &sim->console;

    window_params->window_size_x = (float)event->window.data1;
    window_params->window_size_y = (float)event->window.data2;

    console->cmd_pos_x = 0.02F * window_params->window_size_x;
    console->cmd_pos_y = window_params->window_size_y - (0.1F * window_params->window_size_y);

    console->log_pos_x = console->cmd_pos_x;
    console->log_pos_y = console->cmd_pos_y + (0.05F * window_params->window_size_y);
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// MAIN EVENT CHECKING FUNCTION
////////////////////////////////////////////////////////////////////////////////////////////////////
// the event handling code... checks if events are happening for input and does a task based on that input
void runEventCheck(SDL_Event* event, sim_properties_t* sim, binary_filenames_t* filenames) {
    window_params_t* window_params = &sim->window_params;

    window_params->is_zooming = false; // reset zooming checks
    window_params->is_zooming_in = false;
    window_params->is_zooming_out = false;

    while (SDL_PollEvent(event)) {
        // check if x button is pressed to quit
        if (event->type == SDL_EVENT_QUIT) {
            window_params->reset_sim = true;
            window_params->window_open = false;
            window_params->sim_running = false;
        }
        // check if mouse is moving to update hover state
        else if (event->type == SDL_EVENT_MOUSE_MOTION && event->window.windowID == window_params->main_window_ID) {
            handleMouseMotionEvent(event, sim);
        }
        // check if mouse button is clicked
        else if (event->type == SDL_EVENT_MOUSE_BUTTON_DOWN && event->window.windowID == window_params->main_window_ID) {
            handleMouseButtonDownEvent(event, sim);
        }
        // check if mouse button is released
        else if (event->type == SDL_EVENT_MOUSE_BUTTON_UP && event->window.windowID == window_params->main_window_ID) {
            handleMouseButtonUpEvent(event, sim);
        }
        // check if scroll
        else if (event->type == SDL_EVENT_MOUSE_WHEEL && event->window.windowID == window_params->main_window_ID) {
            handleMouseWheelEvent(event, sim);
        }
        // check if keyboard key is pressed
        else if (event->type == SDL_EVENT_KEY_DOWN && event->window.windowID == window_params->main_window_ID) {
            handleKeyboardEvent(event, sim, filenames);
        }
        // check if text input
        else if (event->type == SDL_EVENT_TEXT_INPUT && event->window.windowID == window_params->main_window_ID) {
            handleTextInputEvent(event, sim);
        }
        // check if window is resized
        else if (event->type == SDL_EVENT_WINDOW_RESIZED &&
                 event->window.windowID == window_params->main_window_ID) {
            handleWindowResizeEvent(event, sim);
        }
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// CONSOLE RENDERING
////////////////////////////////////////////////////////////////////////////////////////////////////
void renderCMDWindow(sim_properties_t sim, font_t* font) {
    console_t* console = &sim.console;
    const window_params_t* window_params = &sim.window_params;

    // display what's in the text box
    addText(font, console->cmd_pos_x, console->cmd_pos_y, console->cmd_text_box, 1.0F);

    // blinking cursor
    if ((window_params->frame_counter / 30) % 2 == 0) {
        char cursor_text[256];
        snprintf(cursor_text, sizeof(cursor_text), "%s_", console->cmd_text_box);
        addText(font, console->cmd_pos_x, console->cmd_pos_y, cursor_text, 1.0F);
    }

    // display whatever text is in the log box (resets when enter is pressed)
    addText(font, console->log_pos_x, console->log_pos_y, console->log, 0.8F);
}
