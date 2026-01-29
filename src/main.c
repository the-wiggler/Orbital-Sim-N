/* //////////////////////////////////////////////////////////////////// *
*      ____  ____  ____  __________     _____ ______  ___     _   __    *
*     / __ \/ __ \/ __ )/  _/_  __/    / ___//  _/  |/  /    / | / /    *
*    / / / / /_/ / __  |/ /  / /       \__ \ / // /|_/ /    /  |/ /     *
*   / /_/ / _, _/ /_/ // /  / /       ___/ // // /  / /    / /|  /      *
*   \____/_/ |_/_____/___/ /_/       /____/___/_/  /_/    /_/ |_/       *
*                                                                       *
*   Author: toastyy-1                                                   *
*   Description: This is the main file for the Orbital Sim N Software.  *
*                                                                       *
* ////////////////////////////////////////////////////////////////////  */


#include <stdio.h>
#include <stdbool.h>

#ifdef _WIN32
    #include <minwindef.h>
#else
    #include <pthread.h>
#endif

#include "globals.h"
#include "types.h"
#include "sim/simulation.h"
#include "utility/telemetry_export.h"
#include "utility/sim_thread.h"

// GUI INCLUDES
#ifdef GUI_ENABLED
#include "gui/SDL_engine.h"
#include "gui/GL_renderer.h"
#include "gui/models.h"

//#include <SDL3/SDL_hints.h> removed because clang tidy got mad... if there's an error on linux try uncommenting this first

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

#include <SDL3/SDL_init.h>
#include <SDL3/SDL_video.h>
#include <SDL3/SDL_events.h>

#endif
// END GUI INCLUDES

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
DWORD WINAPI physicsSim(LPVOID args) {
#else
void* physicsSim(void* args) {
#endif
    sim_properties_t* sim = (sim_properties_t*)args;
    while (sim->window_params.window_open) {
        while (sim->window_params.sim_running) {
            // lock mutex before accessing data
            mutex_lock(&sim_mutex);

            // DOES ALL BODY AND CRAFT CALCULATIONS:
            runCalculations(sim);

            // unlock mutex when done :)
            mutex_unlock(&sim_mutex);
        }
    }
#ifdef _WIN32
    return 0;
#else
    return NULL;
#endif
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// GUI MAIN :)
////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef GUI_ENABLED
int main(int argc, char *argv[]) {
    (void)argc;
    (void)argv;
    ////////////////////////////////////////
    // INIT                               //
    ////////////////////////////////////////
    // initialize simulation objects
    sim_properties_t sim = {
        .global_bodies = {0},
        .global_spacecraft = {0},
        .window_params = {0}
    };

    // CSV file creation
    binary_filenames_t filenames = {
        .global_data_FILE = NULL
    };
#ifdef _WIN32
    fopen_s(&filenames.global_data_FILE, "osn_telem.csv", "w");
#else
    filenames.global_data_FILE = fopen("osn_telem.csv", "w");
#endif
    writeCSVHeader(filenames.global_data_FILE);

#ifdef __linux__
    // force X11 on Linux (fixes SDL text input issues on wayland)
    SDL_SetHint(SDL_HINT_VIDEO_DRIVER, "x11");
#endif

    // initialize SDL
    SDL_Init(SDL_INIT_VIDEO);

    // window parameters & command prompt init
    sim.window_params = init_window_params();
    sim.console = init_console(sim.window_params);

    // SDL and OpenGL window
    SDL_GL_init_t windowInit = init_SDL_OPENGL_window("Orbit Simulation N",
        (int)sim.window_params.window_size_x, (int)sim.window_params.window_size_y, &sim.window_params.main_window_ID);
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
    sphere_mesh_t sphere_mesh;
    sphere_mesh = generateUnitSphere(15, 15);
    VBO_t sphere_buffer = createVBO(sphere_mesh.vertices, sphere_mesh.data_size);
    sim.window_params.planet_model_vertex_count = (int)sphere_mesh.vertex_count; // I couldn't think of a better way to do this ngl

    // create batch to hold all the line geometries we would ever want to draw
    line_batch_t line_batch;
    line_batch = createLineBatch(MAX_LINE_BATCH);

    // craft path tracking
    craft_path_storage_t craft_paths = {0};

    // initialize font for text rendering
    font_t font;
    font = initFont("assets/font.ttf", 24.0F);
    if (font.shader == 0) {
        displayError("Font Error", "Failed to initialize font. Check console for details.");
        return 1;
    }

    ////////////////////////////////////////
    // SIM THREAD INIT                    //
    ////////////////////////////////////////
    // initialize mutex
    mutex_init(&sim_mutex);

#ifdef _WIN32
    #include <winnt.h>
    #include <processthreadsapi.h>
    HANDLE sim_thread = CreateThread(NULL, 0, physicsSim, &sim, 0, NULL);
#else
    // NOLINTNEXTLINE(misc-include-cleaner) - pthread_t from conditional pthread.h include above
    pthread_t simThread;
    pthread_create(&simThread, NULL, physicsSim, &sim);
#endif

    ////////////////////////////////////////////////////////
    // simulation loop                                    //
    ////////////////////////////////////////////////////////
    // default time step
    sim.window_params.time_step = 0.01;

    while (sim.window_params.window_open) {
        // clears previous frame from the screen
        glClearColor(0.0F, 0.0F, 0.0F, 0.0F);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        // user input event checking logic (modifies UI state, no lock needed)
        SDL_Event event;
        runEventCheck(&event, &sim);

        // lock mutex and quickly snapshot simulation data for rendering
        mutex_lock(&sim_mutex);

        // make a quick copy for rendering
        const sim_properties_t sim_copy = sim;

        mutex_unlock(&sim_mutex);

        ////////////////////////////////////////////////////////
        // OPENGL RENDERER
        ////////////////////////////////////////////////////////
        // update viewport for window resizing
        glViewport(0, 0, (int)sim_copy.window_params.window_size_x, (int)sim_copy.window_params.window_size_y);

        // use shader program
        glUseProgram(shaderProgram);

        // casts the camera to the required orientation and zoom (always points to the origin)
        castCamera(sim_copy, shaderProgram);

        // draw coordinate plane
        renderCoordinatePlane(sim_copy, &line_batch);

        // draw planets
        renderPlanets(sim_copy, shaderProgram, sphere_buffer);

        // draw crafts
        renderCrafts(sim_copy, shaderProgram, cone_buffer);

        // stats display
        renderStats(sim_copy, &font);

        // renders visuals things if they are enabled
        renderVisuals(sim_copy, &line_batch, &craft_paths);

        // command window display
        renderCMDWindow(sim_copy, &font);

        // render all queued lines
        renderLines(&line_batch, shaderProgram);

        // render all queued text
        renderText(&font, sim_copy.window_params.window_size_x, sim_copy.window_params.window_size_y, 1, 1, 1);
        ////////////////////////////////////////////////////////
        // END OPENGL RENDERER
        ////////////////////////////////////////////////////////

        // log data
        if (sim.window_params.data_logging_enabled && sim.window_params.sim_running) {
            exportTelemetryCSV(filenames, sim_copy);
        }

        // check if sim needs to be reset
        if (sim.window_params.reset_sim) {
            mutex_lock(&sim_mutex);

            resetSim(&sim);

            mutex_unlock(&sim_mutex);

            // reset paths
            for (int i = 0; i < craft_paths.num_objects; i++) {
                craft_paths.counts[i] = 0;
            }
        }

        // increment frame counter
        sim.window_params.frame_counter++;

        // present the renderer to the screen
        SDL_GL_SwapWindow(window);

#ifdef __EMSCRIPTEN__
        emscripten_sleep(0);
#endif
    }
    ////////////////////////////////////////////////////////
    // end of simulation loop                             //
    ////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////
    // CLEAN UP                                       //
    ////////////////////////////////////////////////////

    // wait for simulation thread
#ifdef _WIN32
    // NOLINTNEXTLINE(misc-include-cleaner)
    WaitForSingleObject(sim_thread, INFINITE);
    // NOLINTNEXTLINE(misc-include-cleaner)
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

    if (filenames.global_data_FILE != NULL) {
        fclose(filenames.global_data_FILE);
    }
    SDL_GL_DestroyContext(glctx);
    SDL_DestroyWindow(window);
    SDL_Quit();
    return 0;
}
#endif


