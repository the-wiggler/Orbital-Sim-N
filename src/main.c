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
#include <GL/gl.h>

#include <stdbool.h>
#include "gui/GL_renderer.h"

// NOTE: ALL CALCULATIONS SHOULD BE DONE IN BASE SI UNITS

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

    // window parameters
    sim.wp = init_window_params();

    // SDL and OpenGL window
    SDL_GL_init_t windowInit = init_SDL_OPENGL_window("Orbit Simulation N",
        (int)sim.wp.window_size_x, (int)sim.wp.window_size_y, &sim.wp.main_window_ID);
    SDL_Window* window = windowInit.window;
    SDL_GLContext glctx = windowInit.glContext;

    // create the shader programs
    GLuint shaderProgram = createShaderProgram("../shaders/simple.vert", "../shaders/simple.frag");

    // enable depth testing
    glEnable(GL_DEPTH_TEST);

    ////////////////////////////////////////
    // CUBE SETUP                         //
    ////////////////////////////////////////
    // create cube vertices with colors (position xyz + color rgb)
    float unit_cube_vertices[] = {
        // positions          // colors
        // Front face
        -0.5f, -0.5f,  0.5f,  1.0f, 0.0f, 0.0f,
         0.5f, -0.5f,  0.5f,  1.0f, 0.0f, 0.0f,
         0.5f,  0.5f,  0.5f,  1.0f, 0.0f, 0.0f,
         0.5f,  0.5f,  0.5f,  1.0f, 0.0f, 0.0f,
        -0.5f,  0.5f,  0.5f,  1.0f, 0.0f, 0.0f,
        -0.5f, -0.5f,  0.5f,  1.0f, 0.0f, 0.0f,

        // Back face
         0.5f, -0.5f, -0.5f,  0.0f, 1.0f, 0.0f,
        -0.5f, -0.5f, -0.5f,  0.0f, 1.0f, 0.0f,
        -0.5f,  0.5f, -0.5f,  0.0f, 1.0f, 0.0f,
        -0.5f,  0.5f, -0.5f,  0.0f, 1.0f, 0.0f,
         0.5f,  0.5f, -0.5f,  0.0f, 1.0f, 0.0f,
         0.5f, -0.5f, -0.5f,  0.0f, 1.0f, 0.0f,

        // Top face
        -0.5f,  0.5f,  0.5f,  0.0f, 0.0f, 1.0f,
         0.5f,  0.5f,  0.5f,  0.0f, 0.0f, 1.0f,
         0.5f,  0.5f, -0.5f,  0.0f, 0.0f, 1.0f,
         0.5f,  0.5f, -0.5f,  0.0f, 0.0f, 1.0f,
        -0.5f,  0.5f, -0.5f,  0.0f, 0.0f, 1.0f,
        -0.5f,  0.5f,  0.5f,  0.0f, 0.0f, 1.0f,

        // Bottom face
        -0.5f, -0.5f, -0.5f,  1.0f, 1.0f, 0.0f,
         0.5f, -0.5f, -0.5f,  1.0f, 1.0f, 0.0f,
         0.5f, -0.5f,  0.5f,  1.0f, 1.0f, 0.0f,
         0.5f, -0.5f,  0.5f,  1.0f, 1.0f, 0.0f,
        -0.5f, -0.5f,  0.5f,  1.0f, 1.0f, 0.0f,
        -0.5f, -0.5f, -0.5f,  1.0f, 1.0f, 0.0f,

        // Right face
         0.5f, -0.5f,  0.5f,  1.0f, 0.0f, 1.0f,
         0.5f, -0.5f, -0.5f,  1.0f, 0.0f, 1.0f,
         0.5f,  0.5f, -0.5f,  1.0f, 0.0f, 1.0f,
         0.5f,  0.5f, -0.5f,  1.0f, 0.0f, 1.0f,
         0.5f,  0.5f,  0.5f,  1.0f, 0.0f, 1.0f,
         0.5f, -0.5f,  0.5f,  1.0f, 0.0f, 1.0f,

        // Left face
        -0.5f, -0.5f, -0.5f,  0.0f, 1.0f, 1.0f,
        -0.5f, -0.5f,  0.5f,  0.0f, 1.0f, 1.0f,
        -0.5f,  0.5f,  0.5f,  0.0f, 1.0f, 1.0f,
        -0.5f,  0.5f,  0.5f,  0.0f, 1.0f, 1.0f,
        -0.5f,  0.5f, -0.5f,  0.0f, 1.0f, 1.0f,
        -0.5f, -0.5f, -0.5f,  0.0f, 1.0f, 1.0f
    };

    float axis_lines[] = {
        // x axis (b
        -10.0f, 0.0f, 0.0f,  0.3f, 0.0f, 0.0f,
         10.0f, 0.0f, 0.0f,  0.3f, 0.0f, 0.0f,

        // y axis (green)
         0.0f, -10.0f, 0.0f,  0.0f, 0.3f, 0.0f,
         0.0f,  10.0f, 0.0f,  0.0f, 0.3f, 0.0f,

        // z axis (blue)
         0.0f, 0.0f, -10.0f,  0.0f, 0.0f, 0.3f,
         0.0f, 0.0f,  10.0f,  0.0f, 0.0f, 0.3f,

        // cross lines
        10.0f, 0.0f, -10.0f,  0.3f, 0.3f, 0.3f,
        -10.0f, 0.0f,  10.0f,  0.3f, 0.3f, 0.3f,

        10.0f, 0.0f, 10.0f,  0.3f, 0.3f, 0.3f,
        -10.0f, 0.0f,  -10.0f,  0.3f, 0.3f, 0.3f,
    };

    VBO_t unit_cube_buffer = createVBO(unit_cube_vertices, sizeof(unit_cube_vertices));

    VBO_t axes_buffer = createVBO(axis_lines, sizeof(axis_lines));

    ////////////////////////////////////////
    // SIM THREAD INIT                    //
    ////////////////////////////////////////
    // initialize simulation thread
    pthread_t simThread;
    pthread_mutex_init(&sim_vars_mutex, NULL);

    // creates the sim thread
    if (pthread_create(&simThread, NULL, physicsSim, &sim) != 0) {
        displayError("ERROR", "Error when creating physics simulation process");
        sim.wp.sim_running = false;
        sim.wp.window_open = false;
        return 1;
    }

    // temp: load JSON for now by default
    readSimulationJSON("simulation_data.json", &sim.gb, &sim.gs);
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

        // use shader progra
        glUseProgram(shaderProgram);

        // casts the camera to the required orientation and zoom (always points to the origin)
        castCamera(sim, shaderProgram);

        // draw coordinate plane
        glBindVertexArray(axes_buffer.VAO);
        mat4 axes_model_mat = mat4_scale(sim.wp.zoom, sim.wp.zoom, sim.wp.zoom);
        setMatrixUniform(shaderProgram, "model", &axes_model_mat);
        glDrawArrays(GL_LINES, 0, 10);

        // draw planets
        glBindVertexArray(unit_cube_buffer.VAO);
        for (int i = 0; i < sim.gb.count; i++) {
            float size_scale_factor = (float)sim.gb.radius[i] / SCALE;
            mat4 scale_mat = mat4_scale(size_scale_factor, size_scale_factor, size_scale_factor);
            mat4 translate_mat = mat4_translation((float)sim.gb.pos_x[i] / SCALE, (float)sim.gb.pos_y[i] / SCALE, (float)sim.gb.pos_z[i] / SCALE);
            mat4 planet_model = mat4_mul(scale_mat, translate_mat);
            setMatrixUniform(shaderProgram, "model", &planet_model);//
            glDrawArrays(GL_TRIANGLES, 0, 36);
        }

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
    glDeleteProgram(shaderProgram);

    fclose(filenames.global_data_FILE);
    SDL_GL_DestroyContext(glctx);
    SDL_DestroyWindow(window);
    SDL_Quit();
    return 0;
}