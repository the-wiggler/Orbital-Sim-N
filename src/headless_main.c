/* //////////////////////////////////////////////////////////////////// *
*      ____  ____  ____  __________     _____ ______  ___     _   __    *
*     / __ \/ __ \/ __ \/  _/_  __/    / ___//  _/  |/  /    / | / /    *
*    / / / / /_/ / __  |/ /  / /       \__ \ / // /|_/ /    /  |/ /     *
*   / /_/ / _, _/ /_/ // /  / /       ___/ // // /  / /    / /|  /      *
*   \____/_/ |_/_____/___/ /_/       /____/___/_/  /_/    /_/ |_/       *
*                                                                       *
*   Author: toastyy-1                                                   *
*   Description: Headless (command-line) main for Orbital Sim N         *
*                                                                       *
* ////////////////////////////////////////////////////////////////////  */

#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "globals.h"
#include "sim/simulation.h"
#include "utility/json_loader.h"
#include "utility/telemetry_export.h"
#include "utility/sim_thread.h"

#ifdef _WIN32
    #include <windows.h>
#else
    #include <pthread.h>
#endif

// Global mutex definition
mutex_t sim_mutex;

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
// HEADLESS ERROR DISPLAY (stub for GUI function)
////////////////////////////////////////////////////////////////////////////////////////////////////
// In headless mode, just print errors to stderr instead of showing a dialog box
void displayError(const char* title, const char* message) {
    fprintf(stderr, "[ERROR] %s: %s\n", title, message);
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// COMMAND LINE INPUT THREAD
////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef _WIN32
DWORD WINAPI commandInputThread(LPVOID args) {
#else
// NOLINTNEXTLINE(readability-function-cognitive-complexity)
void* commandInputThread(void* args) {
#endif
    sim_properties_t* sim = (sim_properties_t*)args;
    char input[256];

    while (sim->window_params.window_open) {
        if (fgets(input, sizeof(input), stdin) != NULL) {
            // remove newline character
            size_t len = strlen(input);
            if (len > 0 && input[len-1] == '\n') {
                input[len-1] = '\0';
            }
            for (size_t i = 0; input[i]; i++) {
                input[i] = tolower(input[i]);
            }

            // parse commands
            if (strcmp(input, "start") == 0 || strcmp(input, "run") == 0) {
                mutex_lock(&sim_mutex);
                sim->window_params.sim_running = true;
                mutex_unlock(&sim_mutex);
                printf("[INFO] Simulation started.\n");
            }
            else if (strcmp(input, "stop") == 0 || strcmp(input, "pause") == 0) {
                mutex_lock(&sim_mutex);
                sim->window_params.sim_running = false;
                mutex_unlock(&sim_mutex);
                printf("[INFO] Simulation paused.\n");
            }
            else if (strcmp(input, "quit") == 0 || strcmp(input, "exit") == 0) {
                mutex_lock(&sim_mutex);
                sim->window_params.sim_running = false;
                sim->window_params.window_open = false;
                mutex_unlock(&sim_mutex);
                printf("[INFO] Shutting down...\n");
                break;
            }
            else if (strcmp(input, "status") == 0) {
                mutex_lock(&sim_mutex);
                printf("[STATUS] Simulation is %s\n", sim->window_params.sim_running ? "RUNNING" : "PAUSED");
                printf("[STATUS] Simulation time: %.2f s\n", sim->window_params.sim_time);
                printf("[STATUS] Time step: %.6f s\n", sim->window_params.time_step);
                printf("[STATUS] Bodies: %d, Spacecraft: %d\n", sim->global_bodies.count, sim->global_spacecraft.count);
                mutex_unlock(&sim_mutex);
            }
            else if (strcmp(input, "reset") == 0) {
                sim->window_params.reset_sim = true;
                printf("[INFO] Simulation reset.\n");
            }
            else if (strcmp(input, "load") == 0) {
                if (sim->global_bodies.count == 0) {
                    printf("[INFO] Loaded sim parameters from JSON\n");
                    readSimulationJSON(SIMULATION_FILENAME, &sim->global_bodies, &sim->global_spacecraft);
                }
                else { printf("[WARNING] Sim JSON already loaded.\n"); }
            }
            else if (strncmp(input, "step ", 5) == 0) {
                // parse the numeric argument after "step "
                char* arg = input + 5;
                // skip leading whitespace
                while (*arg == ' ') { arg++; }

                double new_step = atof(arg);
                if (new_step > 0.0) {
                    mutex_lock(&sim_mutex);
                    sim->window_params.time_step = new_step;
                    mutex_unlock(&sim_mutex);
                    printf("[INFO] Time step changed to %.6f seconds.\n", new_step);
                } else {
                    printf("[ERROR] Invalid time step. Please provide a positive number.\n");
                }
            }
            else if (strncmp(input, "logfreq ", 8) == 0) {
            }
            else if (strcmp(input, "help") == 0) {
                printf("\nAvailable commands:\n");
                printf("  start/run      - Start the simulation\n");
                printf("  stop/pause     - Pause the simulation\n");
                printf("  status         - Show simulation status\n");
                printf("  reset          - Reset the simulation\n");
                printf("  load           - Load data from sim parameter JSON\n");
                printf("  step <num>     - Change the sim time step (e.g., step 0.01)\n");
                printf("  quit/exit      - Exit the program\n");
                printf("  help           - Show this help message\n\n");
            }
            else if (strlen(input) > 0) {
                printf("[ERROR] Unknown command: '%s'. Type 'help' for available commands.\n", input);
            }

            printf("> ");
            fflush(stdout);
        }
    }

#ifdef _WIN32
    return 0;
#else
    return NULL;
#endif
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// HEADLESS MAIN
////////////////////////////////////////////////////////////////////////////////////////////////////
int main(int argc, char *argv[]) {
    (void)argc;
    (void)argv;

    // welcome message
    printf("==================================================================\n\n");
    printf("   ____  ____  ____  __________     _____ ______  ___     _   __\n");
    printf("  / __ \\/ __ \\/ __ )/  _/_  __/    / ___//  _/  |/  /    / | / /\n");
    printf(" / / / / /_/ / __  |/ /  / /       \\__ \\ / // /|_/ /    /  |/ /\n");
    printf("/ /_/ / _, _/ /_/ // /  / /       ___/ // // /  / /    / /|  /\n");
    printf("\\____/_/ |_/_____/___/ /_/       /____/___/_/  /_/    /_/ |_/\n\n");
    printf("==================================================================\n");

    ////////////////////////////////////////
    // INIT                               //
    ////////////////////////////////////////
    // initialize simulation objects
    printf("Initializing simulation objects... ");
    sim_properties_t sim = {
        .global_bodies = {0},
        .global_spacecraft = {0},
        .window_params = {0}
    };
    sim.window_params.window_open = true;  // keep program running
    sim.window_params.sim_running = false; // start paused
    printf("[OK]\n");

    // CSV file creation
    printf("Creating log file... ");
    binary_filenames_t filenames = { .global_data_FILE = fopen("osn_telem.csv", "w") };
    if (!filenames.global_data_FILE) {
        printf("[FAILED]\n");
    } else {
        printf("[OK]\n");
        writeCSVHeader(filenames.global_data_FILE);
    }

    ////////////////////////////////////////
    // SIM THREAD INIT                    //
    ////////////////////////////////////////
    // Initialize mutex
    mutex_init(&sim_mutex);

    printf("Creating simulation process thread... ");
    #ifdef _WIN32
    HANDLE sim_thread = CreateThread(NULL, 0, physicsSim, &sim, 0, NULL);
    if (sim_thread == NULL) {
        printf("[FAILED]\n");
    } else {
        printf("[OK]\n");
    }
    #else
    pthread_t simThread;
    if (pthread_create(&simThread, NULL, physicsSim, &sim)) {
        printf("[FAILED]\n");
    } else {
        printf("[OK]\n");
    }
    #endif

    printf("Creating command input thread... ");
    #ifdef _WIN32
    HANDLE input_thread = CreateThread(NULL, 0, commandInputThread, &sim, 0, NULL);
    if (input_thread == NULL) {
        printf("[FAILED]\n");
    } else {
        printf("[OK]\n");
    }
    #else
    pthread_t inputThread;
    if (pthread_create(&inputThread, NULL, commandInputThread, &sim)) {
        printf("[FAILED]\n");
    } else {
        printf("[OK]\n");
    }
    #endif

    printf("==================================================================\n");
    printf("\nType 'help' for available commands.\n");
    printf("> ");
    fflush(stdout);

    ////////////////////////////////////////////////////////
    // simulation loop                                    //
    ////////////////////////////////////////////////////////
    // default time step
    sim.window_params.time_step = 0.01;
    double csv_update_period = 12960.0F; // updates every n seconds
    double last_csv_update_time = 0;

    while (sim.window_params.window_open) {
        // lock mutex and quickly snapshot simulation data
        mutex_lock(&sim_mutex);

        // make a quick copy for use in main loop
        const sim_properties_t sim_copy = sim;

        mutex_unlock(&sim_mutex);

        // log data on interval -- data logging is enabled by default
        if (sim_copy.window_params.sim_running) {
            double time_since_last_export = sim_copy.window_params.sim_time - last_csv_update_time;
            if (time_since_last_export >= csv_update_period) {
                exportTelemetryCSV(filenames, sim_copy);
                last_csv_update_time = sim_copy.window_params.sim_time;
            }
        }

        // check if sim needs to be reset
        if (sim_copy.window_params.reset_sim) {
            mutex_lock(&sim_mutex);

            resetSim(&sim);

            mutex_unlock(&sim_mutex);
        }
    }
    ////////////////////////////////////////////////////////
    // end of simulation loop                             //
    ////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////
    // CLEAN UP                                       //
    ////////////////////////////////////////////////////

    // wait for threads to complete
    printf("Waiting for threads to finish...\n");
#ifdef _WIN32
    WaitForSingleObject(sim_thread, INFINITE);
    CloseHandle(sim_thread);
    WaitForSingleObject(input_thread, INFINITE);
    CloseHandle(input_thread);
#else
    pthread_join(simThread, NULL);
    pthread_join(inputThread, NULL);
#endif

    // destroy mutex (cross-platform)
    mutex_destroy(&sim_mutex);

    // cleanup all allocated sim memory
    cleanup(&sim);

    if (filenames.global_data_FILE) {
        fclose(filenames.global_data_FILE);
    }

    printf("Goodbye!\n");
    return 0;
}
