#include <stdio.h>
#include <math.h>
#include "globals.h"
#include "types.h"
#include "sim/simulation.h"
#include "gui/SDL_engine.h"
#include "utility/json_loader.h"
#include "utility/telemetry_export.h"
#ifdef _WIN32
    #include <windows.h>
#else
    #include <unistd.h>
#endif
#include <pthread.h>
#include <stdlib.h>
#include <SDL3/SDL.h>
#include <SDL3/SDL_main.h>
#include <GL/glew.h>
#ifdef __APPLE__
#include <OpenGL/gl.h>
#else
#include <GL/gl.h>
#endif
#include <stdbool.h>
#include "gui/GL_renderer.h"
#include "math/matrix.h"
#include "gui/models.h"

// NOTE: ALL CALCULATIONS SHOULD BE DONE IN BASE SI UNITS

// Global mutex definition
pthread_mutex_t sim_vars_mutex;

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

    // enable depth testing
    glEnable(GL_DEPTH_TEST);

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

    ////////////////////////////////////////
    // SIM THREAD INIT                    //
    ////////////////////////////////////////
    // initialize simulation thread
    pthread_t simThread;
    pthread_mutex_init(&sim_vars_mutex, NULL);//

    // creates the sim thread
    if (pthread_create(&simThread, NULL, physicsSim, &sim) != 0) {
        displayError("ERROR", "Error when creating physics simulation process");
        sim.wp.sim_running = false;
        sim.wp.window_open = false;
        return 1;
    }

    // default time step
    sim.wp.time_step = 0.01;

    ////////////////////////////////////////////////////////
    // simulation loop                                    //
    ////////////////////////////////////////////////////////
    while (sim.wp.window_open) {
        // clears previous frame from the screen
        glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        // lock body_sim mutex for reading
        pthread_mutex_lock(&sim_vars_mutex);

        // user input event checking logic
        SDL_Event event;
        runEventCheck(&event, &sim);

        ////////////////////////////////////////////////////////
        // OPENGL RENDERER
        ////////////////////////////////////////////////////////
        // update viewport for window resizing
        glViewport(0, 0, (int)sim.wp.window_size_x, (int)sim.wp.window_size_y);

        // use shader program
        glUseProgram(shaderProgram);

        // casts the camera to the required orientation and zoom (always points to the origin)
        castCamera(sim, shaderProgram);

        // draw coordinate plane
        renderCoordinatePlane(sim, &line_batch);

        // draw planets
        renderPlanets(sim, shaderProgram, sphere_buffer);

        // draw crafts
        renderCrafts(sim, shaderProgram, cone_buffer);

        // stats display
        renderStats(sim, &font);

        // renders debug things if they are enabled
        renderDebug(sim, &line_batch);

        // command window display
        renderCMDWindow(&sim, &font);

        // render all queued lines
        renderLines(&line_batch, shaderProgram);

        // render all queued text
        renderText(&font, sim.wp.window_size_x, sim.wp.window_size_y, 1, 1, 1);
        ////////////////////////////////////////////////////////
        // END OPENGL RENDERER
        ////////////////////////////////////////////////////////

        if (sim.wp.data_logging_enabled) {
            exportTelemetryBinary(filenames, &sim);
        }

        // check if sim needs to be reset
        if (sim.wp.reset_sim) resetSim(&sim);

        // unlock sim vars mutex when done
        pthread_mutex_unlock(&sim_vars_mutex);

        // increment frame counter
        sim.wp.frame_counter++;

        // present the renderer to the screen
        SDL_GL_SwapWindow(window);
    }
    ////////////////////////////////////////////////////////
    // end of simulation loop                             //
    ////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////
    // CLEAN UP                                       //
    ////////////////////////////////////////////////////

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