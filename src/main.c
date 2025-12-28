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
#include "math/matrix.h"
#include "gui/models.h"

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

    // initialize the text renderer
    text_renderer_t textRenderer = initTextRenderer("font.ttf", 24, (int)sim.wp.window_size_x, (int)sim.wp.window_size_y);
    float white_color[] = { 1.0f, 1.0f, 1.0f};

    ////////////////////////////////////////
    // MESH/BUFFER SETUP                  //
    ////////////////////////////////////////

    VBO_t unit_cube_buffer = createVBO(UNIT_CUBE_VERTICES, sizeof(UNIT_CUBE_VERTICES));

    VBO_t axes_buffer = createVBO(AXES_LINES, sizeof(AXES_LINES));

    VBO_t cone_buffer = createVBO(CONE_VERTICES, sizeof(CONE_VERTICES));

    sphere_mesh_t sphere_mesh = generateUnitSphere(15, 15);
    VBO_t sphere_buffer = createVBO(sphere_mesh.vertices, sphere_mesh.data_size);

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

        // use shader program
        glUseProgram(shaderProgram);

        // casts the camera to the required orientation and zoom (always points to the origin)
        castCamera(sim, shaderProgram);

        // draw coordinate plane
        glBindVertexArray(axes_buffer.VAO);
        mat4 axes_model_mat = mat4_scale(sim.wp.zoom, sim.wp.zoom, sim.wp.zoom);
        setMatrixUniform(shaderProgram, "model", &axes_model_mat);
        glDrawArrays(GL_LINES, 0, 10);

        // draw planets
        glBindVertexArray(sphere_buffer.VAO);
        for (int i = 0; i < sim.gb.count; i++) {
            // create a scale matrix based on the radius of the planet (pulled from JSON data)
            float size_scale_factor = (float)sim.gb.radius[i] / SCALE;
            mat4 scale_mat = mat4_scale(size_scale_factor, size_scale_factor, size_scale_factor);
            // create a translation matrix based on the current in-sim-world position of the planet
            mat4 translate_mat = mat4_translation((float)sim.gb.pos_x[i] / SCALE, (float)sim.gb.pos_y[i] / SCALE, (float)sim.gb.pos_z[i] / SCALE);
            // multiply the two matrices together to get a final scale/position matrix for the planet model on the screen
            mat4 planet_model = mat4_mul(translate_mat, scale_mat);
            // apply matrix and render to screen
            setMatrixUniform(shaderProgram, "model", &planet_model);
            glDrawArrays(GL_TRIANGLES, 0, (GLsizei)sphere_mesh.vertex_count);
        }

        // draw crafts
        glBindVertexArray(cone_buffer.VAO);
        for (int i = 0; i < sim.gs.count; i++) {
            // create a scale matrix
            float size_scale_factor = 0.1f; // arbitrary scale based on what looks nice on screen
            mat4 scale_mat = mat4_scale(size_scale_factor, size_scale_factor * 2, size_scale_factor);
            // create a translation matrix based on the current in-sim-world position of the spacecraft
            mat4 translate_mat = mat4_translation((float)sim.gs.pos_x[i] / SCALE, (float)sim.gs.pos_y[i] / SCALE, (float)sim.gs.pos_z[i] / SCALE);
            // multiply the two matrices together to get a final scale/position matrix for the planet model on the screen
            mat4 spacecraft_model = mat4_mul(translate_mat, scale_mat);
            // apply matrix and render to screen
            setMatrixUniform(shaderProgram, "model", &spacecraft_model);
            glDrawArrays(GL_TRIANGLES, 0, 48);
        }
        
        renderText(&textRenderer, "Hello, world!", 10.0f, 10.0f, 1.0f, white_color);

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
    freeSphere(&sphere_mesh);
    deleteVBO(unit_cube_buffer);
    deleteVBO(axes_buffer);
    deleteVBO(sphere_buffer);
    glDeleteProgram(shaderProgram);

    fclose(filenames.global_data_FILE);
    SDL_GL_DestroyContext(glctx);
    SDL_DestroyWindow(window);
    SDL_Quit();
    return 0;
}