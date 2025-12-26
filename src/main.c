#include <stdio.h>
#include <math.h>

#include "globals.h"
#include "types.h"
#include "sim/simulation.h"
#include "gui/SDL_engine.h"
#include "gui/craft_view.h"
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
    ///
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

    // openGL init
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MAJOR_VERSION, 3);
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MINOR_VERSION, 3);
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_PROFILE_MASK, SDL_GL_CONTEXT_PROFILE_CORE);                // core profile
    SDL_GL_SetAttribute(SDL_GL_DOUBLEBUFFER, 1);
    SDL_GL_SetAttribute(SDL_GL_DEPTH_SIZE, 24);
    SDL_GL_SetAttribute(SDL_GL_STENCIL_SIZE, 8);

    // initialize window parameters
    SDL_Init(SDL_INIT_VIDEO);
    init_window_params(&sim.wp);

    // initialize SDL3 window
    // create an SDL window
    SDL_Window* window = SDL_CreateWindow("Orbit Simulation",
        (int)sim.wp.window_size_x, (int)sim.wp.window_size_y,
        SDL_WINDOW_OPENGL | SDL_WINDOW_RESIZABLE);

    // store the window ID for event handling
    sim.wp.main_window_ID = SDL_GetWindowID(window);

    // initialize OpenGL and GLEW
    SDL_GLContext glctx = SDL_GL_CreateContext(window);
    glewExperimental = GL_TRUE;
    glewInit();

    // VSync
    SDL_GL_SetSwapInterval(1);

    printf("OpenGL version: %p\n", glGetString(GL_VERSION));
    printf("GLEW version: %p\n", glewGetString(GLEW_VERSION));

    ////////////////////////////////////////
    // SHADER SETUP                       //
    ////////////////////////////////////////
    // load shader sources from files
    char* vertexShaderSource = loadShaderSource("../shaders/simple.vert");
    char* fragmentShaderSource = loadShaderSource("../shaders/simple.frag");

    // compile the shaders
    GLuint vertexShader = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vertexShader, 1, (const char**)&vertexShaderSource, NULL);
    glCompileShader(vertexShader);

    GLuint fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fragmentShader, 1, (const char**)&fragmentShaderSource, NULL);
    glCompileShader(fragmentShader);

    // link the shaders into a program
    GLuint shaderProgram = glCreateProgram();
    glAttachShader(shaderProgram, vertexShader);
    glAttachShader(shaderProgram, fragmentShader);
    glLinkProgram(shaderProgram);
    glDeleteShader(vertexShader);
    glDeleteShader(fragmentShader);

    // free the shader source strings (bye bye!)
    free(vertexShaderSource);
    free(fragmentShaderSource);

    // circle buffer for rendering
    // this is updated each frame with transformed coordinates
    GLuint VAO_circle, VBO_circle;
    glGenVertexArrays(1, &VAO_circle);
    glGenBuffers(1, &VBO_circle);
    glBindVertexArray(VAO_circle);
    glBindBuffer(GL_ARRAY_BUFFER, VBO_circle);

    // allocate buffer for circle vertices (center + perimeter points)
    #define CIRCLE_SEGMENTS 32
    glBufferData(GL_ARRAY_BUFFER, sizeof(float) * 2 * (CIRCLE_SEGMENTS + 2), NULL, GL_DYNAMIC_DRAW);
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 2 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);

    glBindVertexArray(0);

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
    sim.wp.time_step = 0.0001;

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

        float ndc_x1, ndc_y1, ndc_x2, ndc_y2;

        // transform first planet point in world space
        worldToNDC(sim.gb.pos_x[0], sim.gb.pos_y[0],
                   sim.wp.screen_origin_x, sim.wp.screen_origin_y,
                   sim.wp.meters_per_pixel,
                   sim.wp.window_size_x, sim.wp.window_size_y,
                   &ndc_x1, &ndc_y1);

        // transform second planet point in world space
        worldToNDC(sim.gb.pos_x[1], sim.gb.pos_y[1],
                   sim.wp.screen_origin_x, sim.wp.screen_origin_y,
                   sim.wp.meters_per_pixel,
                   sim.wp.window_size_x, sim.wp.window_size_y,
                   &ndc_x2, &ndc_y2);

        glUseProgram(shaderProgram);
        GLint colorLoc = glGetUniformLocation(shaderProgram, "color");

        // draw circle at first planet
        float circle_vertices[2 * (CIRCLE_SEGMENTS + 2)];

        // Convert world-space radius to NDC coordinates (scales with zoom)
        double world_radius = 1e7;  // radius in meters
        double pixel_radius = world_radius / sim.wp.meters_per_pixel;
        float circle_radius = (float)(2.0 * pixel_radius / sim.wp.window_size_x);

        // center vertex
        circle_vertices[0] = ndc_x1;
        circle_vertices[1] = ndc_y1;

        // perimeter vertices
        for (int i = 0; i <= CIRCLE_SEGMENTS; i++) {
            float angle = 2.0f * 3.14159265f * i / CIRCLE_SEGMENTS;
            circle_vertices[2 + i * 2] = ndc_x1 + circle_radius * cosf(angle);
            circle_vertices[3 + i * 2] = ndc_y1 + circle_radius * sinf(angle);
        }

        glBindBuffer(GL_ARRAY_BUFFER, VBO_circle);
        glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(circle_vertices), circle_vertices);
        glUniform4f(colorLoc, 0.0f, 1.0f, 0.0f, 1.0f); // Green
        glBindVertexArray(VAO_circle);
        glDrawArrays(GL_TRIANGLE_FAN, 0, CIRCLE_SEGMENTS + 2);

        // draw circle at second planet
        circle_vertices[0] = ndc_x2;
        circle_vertices[1] = ndc_y2;

        for (int i = 0; i <= CIRCLE_SEGMENTS; i++) {
            float angle = 2.0f * 3.14159265f * i / CIRCLE_SEGMENTS;
            circle_vertices[2 + i * 2] = ndc_x2 + circle_radius * cosf(angle);
            circle_vertices[3 + i * 2] = ndc_y2 + circle_radius * sinf(angle);
        }

        glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(circle_vertices), circle_vertices);
        glDrawArrays(GL_TRIANGLE_FAN, 0, CIRCLE_SEGMENTS + 2);

        glBindVertexArray(0);


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

    // wait for simulation thread to finishames.global_data_FILE);
    pthread_join(simThread, NULL);

    // destroy mutex
    pthread_mutex_destroy(&sim_vars_mutex);

    // cleanup all allocated memory
    cleanup(&sim);

    // Cleanup OpenGL resources
    glDeleteVertexArrays(1, &VAO_circle);
    glDeleteBuffers(1, &VBO_circle);
    glDeleteProgram(shaderProgram);

    fclose(filenames.global_data_FILE);
    SDL_GL_DestroyContext(glctx);
    SDL_DestroyWindow(window);
    SDL_Quit();
    return 0;
}