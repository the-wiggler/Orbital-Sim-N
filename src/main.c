/* //////////////////////////////////////////////////////////////////// *
*      ____  ____  ____  __________     _____ ______  ___     _   __    *
*     / __ \/ __ \/ __ )/  _/_  __/    / ___//  _/  |/  /    / | / /    *
*    / / / / /_/ / __  |/ /  / /       \__ \ / // /|_/ /    /  |/ /     *
*   / /_/ / _, _/ /_/ // /  / /       ___/ // // /  / /    / /|  /      *
*   \____/_/ |_/_____/___/ /_/       /____/___/_/  /_/    /_/ |_/       *
*                                                                       *
*   Author: toastyy-1                                                   *
*   This is the main file for the Orbital Sim N Software.               *
*                                                                       *
* ////////////////////////////////////////////////////////////////////  */


#include <stdio.h>

#ifdef _WIN32
    #include <minwindef.h>
#else
    #include <pthread.h>
#endif

#ifdef GUI_ENABLED
#ifdef __EMSCRIPTEN__
    #include <emscripten/emscripten.h>
    #include <GLES3/gl3.h>
#else
    #include <GL/glew.h>
#endif

#ifdef __APPLE__
    #include <OpenGL/gl.h>
#else
    #include <GL/gl.h>
#endif

#include <SDL3/SDL_init.h>
#include <SDL3/SDL_video.h>
#include <SDL3/SDL_events.h>
#endif

#include "types.h"
#include "sim/simulation.h"
#include "utility/telemetry_export.h"
#include "utility/sim_thread.h"

#ifdef GUI_ENABLED
#include "gui/SDL_engine.h"
#include "gui/GL_renderer.h"
#endif
// END GUI INCLUDES

// Global mutex definition
mutex_t sim_mutex;
// this is purposely made a global var in this file
// as it is expected that mutex locks should not be
// hidden within other files

// NOTE: ALL CALCULATIONS SHOULD BE DONE IN BASE SI UNITS

////////////////////////////////////////////////////////////////////////////////////////////////////
// PHYSICS SIMULATION THREAD!
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
int main() {
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
        .csv_update_period = 60.0F,
        .last_csv_update_time = 0,
        .global_data_FILE = NULL
    };
#ifdef _WIN32
    fopen_s(&filenames.global_data_FILE, "osn_telem.csv", "w");
#else
    filenames.global_data_FILE = fopen("osn_telem.csv", "w");
#endif
    writeCSVHeader(filenames.global_data_FILE);

    // SDL and OpenGL window
    SDL_GL_init_t windowInit = init_SDL_OPENGL_window("Orbit Simulation N", &sim);
    SDL_Window* window = windowInit.window;
    SDL_GLContext glctx = windowInit.glContext;

    GLuint shaderProgram = init_GL_shader();

    GL_assets_t assets = init_GL_assets(&sim);

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

        // lock mutex and quickly snapshot simulation data for rendering
        mutex_lock(&sim_mutex);

        // user input event checking logic
        SDL_Event event;
        runEventCheck(&event, &sim, &filenames);

        // updates the orbital elements for each orbital body and craft (for visual guidelines)
        updateSystemOrbitalElements(&sim);

        // calculate the total energy of the system and store it for external use
        sim.global_bodies.total_system_energy = calculateTotalSystemEnergy(&sim);

        // update camera target if tracking a spacecraft
        if (sim.window_params.track_craft_id >= 0 && sim.window_params.track_craft_id < sim.global_spacecraft.count) {
            const vec3 craft_pos = sim.global_spacecraft.spacecraft[sim.window_params.track_craft_id].pos;
            sim.window_params.cam_target.x = (float)(craft_pos.x / SCALE);
            sim.window_params.cam_target.y = (float)(craft_pos.y / SCALE);
            sim.window_params.cam_target.z = (float)(craft_pos.z / SCALE);
        }

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
        renderCoordinatePlane(sim_copy, &assets.lines);

        // draw planets
        renderPlanets(sim_copy, shaderProgram, assets.sphere);

        // draw crafts
        renderCrafts(sim_copy, shaderProgram, assets.cone);

        // stats display
        renderStats(sim_copy, &assets.font);

        // renders visuals things if they are enabled
        renderVisuals(sim_copy, &assets.lines, &assets.craft_paths);

        // command window display
        renderCMDWindow(sim_copy, &assets.font);

        // render all queued lines
        renderLines(&assets.lines, shaderProgram);

        // render all queued text
        renderText(&assets.font, sim_copy.window_params.window_size_x, sim_copy.window_params.window_size_y, 1, 1, 1);
        ////////////////////////////////////////////////////////
        // END OPENGL RENDERER
        ////////////////////////////////////////////////////////

        // log data on interval
        if (sim_copy.window_params.data_logging_enabled && sim_copy.window_params.sim_running) {
            double time_since_last_export = sim_copy.window_params.sim_time - filenames.last_csv_update_time;
            if (time_since_last_export >= filenames.csv_update_period) {
                exportTelemetryCSV(filenames, sim_copy);
                filenames.last_csv_update_time = sim_copy.window_params.sim_time;
            }
        }

        // check if sim needs to be reset
        if (sim_copy.window_params.reset_sim) {
            mutex_lock(&sim_mutex);

            resetSim(&sim);

            mutex_unlock(&sim_mutex);

            // reset paths
            for (int i = 0; i < assets.craft_paths.num_objects; i++) {
                assets.craft_paths.counts[i] = 0;
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

    // cleanup OpenGL resources
    glDeleteProgram(shaderProgram);

    SDL_GL_DestroyContext(glctx);
    SDL_DestroyWindow(window);
    SDL_Quit();
    return 0;
}
#endif


