#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <pthread.h>

#include "globals.h"
#include "types.h"
#include "sim/simulation.h"
#include "gui/SDL_engine.h"
#include "utility/telemetry_export.h"
#include "gui/GL_renderer.h"
#include "math/matrix.h"
#include "gui/models.h"

#ifdef _WIN32
    #include <windows.h>
#endif

#include <SDL3/SDL.h>
#include <SDL3/SDL_main.h>

#ifdef __APPLE__
    #include <OpenGL/gl.h>
#else
    #include <GL/gl.h>
#endif

#ifdef __EMSCRIPTEN__
    #include <emscripten/emscripten.h>
    #include <GLES3/gl3.h>
#else
    #include <GL/glew.h>
#endif

// NOTE: ALL CALCULATIONS SHOULD BE DONE IN BASE SI UNITS

// Global mutex definition
pthread_mutex_t sim_vars_mutex;

// Context structure to pass data to the main loop function
// This unifies the requirements for both Native and WASM
typedef struct {
    sim_properties_t* sim;
    SDL_Window* window;
    GLuint shaderProgram;
    VBO_t sphere_buffer;
    VBO_t cone_buffer;
    VBO_t unit_cube_buffer; // currently unused but kept for parity
    line_batch_t* line_batch;
    font_t* font;
    binary_filenames_t* filenames;
} AppContext;

////////////////////////////////////////////////////////////////////////////////////////////////////
// PHYSICS SIMULATION THREAD
////////////////////////////////////////////////////////////////////////////////////////////////////
void* physicsSim(void* args) {
    sim_properties_t* sim = (sim_properties_t*)args;
    while (sim->wp.window_open) {
        while (sim->wp.sim_running) {
            // lock mutex before accessing data
            pthread_mutex_lock(&sim_vars_mutex);

            // DOES ALL BODY AND CRAFT CALCULATIONS:
            runCalculations(sim);

            // unlock mutex when done :)
            pthread_mutex_unlock(&sim_vars_mutex);
        }
    }
    return NULL;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// UNIFIED MAIN LOOP STEP
// This function runs exactly one frame of rendering and logic.
////////////////////////////////////////////////////////////////////////////////////////////////////
void mainLoopStep(void* arg) {
    AppContext* ctx = (AppContext*)arg;
    sim_properties_t* s = ctx->sim;

    // handle loop cancellation for Emscripten if window is closed
    if (!s->wp.window_open) {
        #ifdef __EMSCRIPTEN__
        emscripten_cancel_main_loop();
        #endif
        return;
    }

    // clears previous frame from the screen
    glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // lock body_sim mutex for reading
    pthread_mutex_lock(&sim_vars_mutex);

    // user input event checking logic
    SDL_Event event;
    runEventCheck(&event, s);

    ////////////////////////////////////////////////////////
    // OPENGL RENDERER
    ////////////////////////////////////////////////////////
    // update viewport for window resizing
    glViewport(0, 0, (int)s->wp.window_size_x, (int)s->wp.window_size_y);

    // use shader program
    glUseProgram(ctx->shaderProgram);

    // casts the camera to the required orientation and zoom (always points to the origin)
    castCamera(*s, ctx->shaderProgram);

    // draw coordinate plane
    renderCoordinatePlane(*s, ctx->line_batch);

    // draw planets
    renderPlanets(*s, ctx->shaderProgram, ctx->sphere_buffer);

    // draw crafts
    renderCrafts(*s, ctx->shaderProgram, ctx->cone_buffer);

    // stats display
    renderStats(*s, ctx->font);

    // renders visuals things if they are enabled
    renderVisuals(s, ctx->line_batch);

    // command window display
    renderCMDWindow(s, ctx->font);

    // render all queued lines
    renderLines(ctx->line_batch, ctx->shaderProgram);

    // render all queued text
    renderText(ctx->font, s->wp.window_size_x, s->wp.window_size_y, 1, 1, 1);
    ////////////////////////////////////////////////////////
    // END OPENGL RENDERER
    ////////////////////////////////////////////////////////

    if (s->wp.data_logging_enabled) {
        exportTelemetryBinary(*ctx->filenames, s);
    }

    // check if sim needs to be reset
    if (s->wp.reset_sim) resetSim(s);

    // unlock sim vars mutex when done
    pthread_mutex_unlock(&sim_vars_mutex);

    // increment frame counter
    s->wp.frame_counter++;

    // present the renderer to the screen
    SDL_GL_SwapWindow(ctx->window);
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// MAIN :)
////////////////////////////////////////////////////////////////////////////////////////////////////
int main(int argc, char *argv[]) {
    (void)argc;
    (void)argv;
    ////////////////////////////////////////
    // INIT                               //
    ////////////////////////////////////////
    // initialize simulation objects
    sim_properties_t sim = {
        .gb = {0},
        .gs = {0},
        .wp = {0}
    };

    // binary file creation
    binary_filenames_t filenames = {
        .global_data_FILE = fopen("global_data.bin", "wb")
    };

    // force X11 on Linux (fixes text input issues on wayland)
#ifdef __linux__
    SDL_SetHint(SDL_HINT_VIDEO_DRIVER, "x11");
#endif

    // initialize SDL
    SDL_Init(SDL_INIT_VIDEO);

    // window parameters & command prompt init
    sim.wp = init_window_params();
    sim.console = init_console(sim.wp);

    // SDL and OpenGL window
    SDL_GL_init_t windowInit = init_SDL_OPENGL_window("Orbit Simulation N",
        (int)sim.wp.window_size_x, (int)sim.wp.window_size_y, &sim.wp.main_window_ID);
    SDL_Window* window = windowInit.window;
    SDL_GLContext glctx = windowInit.glContext;

    // create the shader programs
    GLuint shaderProgram = createShaderProgram("shaders/simple.vert", "shaders/simple.frag");
    if (shaderProgram == 0) {
        displayError("Shader Error", "Failed to create shader program. Check console for details.");
        return 1;
    }

    // enable depth testing
    glEnable(GL_DEPTH_TEST);

    // enable blending for transparency
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    ////////////////////////////////////////
    // MESH/BUFFER SETUP                  //
    ////////////////////////////////////////

    // create buffer for cube shape
    VBO_t unit_cube_buffer = createVBO(UNIT_CUBE_VERTICES, sizeof(UNIT_CUBE_VERTICES));

    // create buffer for cone shape
    VBO_t cone_buffer = createVBO(CONE_VERTICES, sizeof(CONE_VERTICES));

    // create buffer for sphere shape
    sphere_mesh_t sphere_mesh = generateUnitSphere(15, 15);
    VBO_t sphere_buffer = createVBO(sphere_mesh.vertices, sphere_mesh.data_size);
    sim.wp.planet_model_vertex_count = (int)sphere_mesh.vertex_count; // I couldn't think of a better way to do this ngl

    // create batch to hold all the line geometries we would ever want to draw!
    line_batch_t line_batch = createLineBatch(1000);

    // initialize font for text rendering
    font_t font = initFont("assets/font.ttf", 24.0f);
    if (font.shader == 0) {
        displayError("Font Error", "Failed to initialize font. Check console for details.");
        return 1;
    }

    ////////////////////////////////////////
    // SIM THREAD INIT                    //
    ////////////////////////////////////////
    // initialize simulation thread
    pthread_t simThread;
    pthread_mutex_init(&sim_vars_mutex, NULL);

    // creates the sim thread
    if (pthread_create(&simThread, NULL, physicsSim, &sim) != 0) {
        displayError("ERROR", "Error when creating simulation process thread");
        sim.wp.sim_running = false;
        sim.wp.window_open = false;
        return 1;
    }

    // default time step
    sim.wp.time_step = 0.01;

    ////////////////////////////////////////////////////////
    // PREPARE APP CONTEXT                                //
    ////////////////////////////////////////////////////////
    AppContext ctx = {
        .sim = &sim,
        .window = window,
        .shaderProgram = shaderProgram,
        .sphere_buffer = sphere_buffer,
        .cone_buffer = cone_buffer,
        .unit_cube_buffer = unit_cube_buffer,
        .line_batch = &line_batch,
        .font = &font,
        .filenames = &filenames
    };

    ////////////////////////////////////////////////////////
    // SIMULATION LOOP                                    //
    ////////////////////////////////////////////////////////

#ifdef __EMSCRIPTEN__
    // WebAssembly: Hand control to the browser's main loop
    emscripten_set_main_loop_arg(mainLoopStep, &ctx, 0, 1);
#else
    // Native: Run standard while loop
    while (sim.wp.window_open) {
        mainLoopStep(&ctx);
    }
#endif

    ////////////////////////////////////////////////////
    // CLEAN UP                                       //
    ////////////////////////////////////////////////////
    // Note: Emscripten usually does not reach here unless main loop is cancelled,
    // but the OS/Browser reclaims memory anyway.

    // wait for simulation thread
    pthread_join(simThread, NULL);

    // destroy mutex
    pthread_mutex_destroy(&sim_vars_mutex);

    // cleanup all allocated sim memory
    cleanup(&sim);

    // Cleanup OpenGL resources
    freeSphere(&sphere_mesh);
    deleteVBO(unit_cube_buffer);
    deleteVBO(sphere_buffer);
    freeLines(&line_batch);
    freeFont(&font);
    glDeleteProgram(shaderProgram);

    fclose(filenames.global_data_FILE);
    SDL_GL_DestroyContext(glctx);
    SDL_DestroyWindow(window);
    SDL_Quit();
    return 0;
}