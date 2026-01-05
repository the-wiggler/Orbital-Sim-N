#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>

#include "globals.h"
#include "types.h"
#include "math/matrix.h"
#include "sim/simulation.h"
#include "gui/SDL_engine.h"
#include "gui/GL_renderer.h"
#include "gui/models.h"
#include "utility/telemetry_export.h"
#include "utility/sim_thread.h"

#ifdef _WIN32
    #include <windows.h>
#else
    #include <pthread.h>
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
mutex_t sim_mutex;
// this is purposely made a global var in this file
// as it is expected that mutex locks should not be
// hidden within other files

// NOTE: ALL CALCULATIONS SHOULD BE DONE IN BASE SI UNITS

////////////////////////////////////////////////////////////////////////////////////////////////////
// PHYSICS SIMULATION THREAD
////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef _WIN32
DWORD WINAPI physicsSim_win32(LPVOID args) {
    sim_properties_t* sim = (sim_properties_t*)args;
    while (sim->wp.window_open) {
        while (sim->wp.sim_running) {
            // lock mutex before accessing data
            mutex_lock(&sim_mutex);

            // DOES ALL BODY AND CRAFT CALCULATIONS:
            runCalculations(sim);

            // unlock mutex when done :)
            mutex_unlock(&sim_mutex);
        }
    }
    return 0;
}
#else
void* physicsSim(void* args) {
    sim_properties_t* sim = (sim_properties_t*)args;
    while (sim->wp.window_open) {
        while (sim->wp.sim_running) {
            // lock mutex before accessing data
            mutex_lock(&sim_mutex);

            // DOES ALL BODY AND CRAFT CALCULATIONS:
            runCalculations(sim);

            // unlock mutex when done :)
            mutex_unlock(&sim_mutex);
        }
    }
    return NULL;
}
#endif

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

#ifdef __linux__
    // force X11 on Linux (fixes text input issues on wayland)
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
    // Initialize mutex
    mutex_init(&sim_mutex);

#ifdef _WIN32
    HANDLE sim_thread = CreateThread(NULL, 0, physicsSim_win32, &sim, 0, NULL);
#else
    pthread_t simThread;
    pthread_create(&simThread, NULL, physicsSim, &sim);
#endif

    ////////////////////////////////////////////////////////
    // SIMULATION LOOP                                    //
    ////////////////////////////////////////////////////////
    // default time step
    sim.wp.time_step = 0.01;

    // Native: Run standard while loop
    while (sim.wp.window_open) {
        // clears previous frame from the screen
        glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        // user input event checking logic (modifies UI state, no lock needed)
        SDL_Event event;
        runEventCheck(&event, &sim);
        // lock mutex and quickly snapshot simulation data for rendering
        mutex_lock(&sim_mutex);
        // make a quick copy for rendering
        sim_properties_t sim_snapshot = sim;
        mutex_unlock(&sim_mutex);
        ////////////////////////////////////////////////////////
        // OPENGL RENDERER
        ////////////////////////////////////////////////////////
        // update viewport for window resizing
        glViewport(0, 0, (int)sim_snapshot.wp.window_size_x, (int)sim_snapshot.wp.window_size_y);
        // use shader program
        glUseProgram(shaderProgram);
        // casts the camera to the required orientation and zoom (always points to the origin)
        castCamera(sim_snapshot, shaderProgram);
        // draw coordinate plane
        renderCoordinatePlane(sim_snapshot, &line_batch);
        // draw planets
        renderPlanets(sim_snapshot, shaderProgram, sphere_buffer);
        // draw crafts
        renderCrafts(sim_snapshot, shaderProgram, cone_buffer);
        // stats display
        renderStats(sim_snapshot, &font);
        // renders visuals things if they are enabled
        renderVisuals(&sim_snapshot, &line_batch);
        // command window display
        renderCMDWindow(&sim_snapshot, &font);
        // render all queued lines
        renderLines(&line_batch, shaderProgram);
        // render all queued text
        renderText(&font, sim_snapshot.wp.window_size_x, sim_snapshot.wp.window_size_y, 1, 1, 1);
        // log data
        if (sim.wp.data_logging_enabled) {
            mutex_lock(&sim_mutex);

            exportTelemetryBinary(filenames, &sim);

            mutex_unlock(&sim_mutex);
        }
        // check if sim needs to be reset
        if (sim.wp.reset_sim) {
            mutex_lock(&sim_mutex);

            resetSim(&sim);

            mutex_unlock(&sim_mutex);
        }
        // increment frame counter
        sim.wp.frame_counter++;
        // present the renderer to the screen
        SDL_GL_SwapWindow(window);

#ifdef __EMSCRIPTEN__
        emscripten_sleep(0);
#endif
    }

    ////////////////////////////////////////////////////
    // CLEAN UP                                       //
    ////////////////////////////////////////////////////
    // Note: Emscripten usually does not reach here unless main loop is cancelled,
    // but the OS/Browser reclaims memory anyway.

    // wait for simulation thread
#ifdef _WIN32
    WaitForSingleObject(sim_thread, INFINITE);
    CloseHandle(sim_thread);
#else
    pthread_join(simThread, NULL);
#endif

    // destroy mutex (cross-platform)
    mutex_destroy(&sim_mutex);

    // cleanup all allocated sim memory
    cleanup(&sim);

    // cleanup OpenGL resources
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