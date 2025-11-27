#ifndef SIM_CALCULATIONS_H
#define SIM_CALCULATIONS_H

#include <SDL3_ttf/SDL_ttf.h>
#include "config.h"

extern const double G;
extern TTF_Font* g_font;

// orbital body stuff
void body_calculateGravForce(body_properties_t* b, body_properties_t b2);
void body_updateMotion(body_properties_t* b, double dt);
void body_transformCoordinates(body_properties_t* b, window_params_t window_params);
void body_calculateKineticEnergy(body_properties_t* b);
float body_calculateVisualRadius(body_properties_t* body, window_params_t wp);
void body_addOrbitalBody(body_properties_t** gb, int* num_bodies, char* name, double mass, double x_pos, double y_pos, double x_vel, double y_vel);
void resetSim(double* sim_time, body_properties_t** gb, int* num_bodies, spacecraft_properties_t** sc, int* num_craft);

// system energy calculation
double calculateTotalSystemEnergy(body_properties_t* gb, spacecraft_properties_t* sc, int num_bodies, int num_craft);

// spacecraft body stuff
void craft_calculateGravForce(spacecraft_properties_t* s, body_properties_t b);
void craft_transformCoordinates(spacecraft_properties_t* s, window_params_t wp);
void craft_updateMotion(spacecraft_properties_t* s, double dt);
void craft_applyThrust(spacecraft_properties_t* s);
void craft_consumeFuel(spacecraft_properties_t* s, double dt);
void craft_addSpacecraft(spacecraft_properties_t** sc, int* num_craft, char* name,
                        double x_pos, double y_pos, double x_vel, double y_vel,
                        double dry_mass, double fuel_mass, double thrust,
                        double specific_impulse, double mass_flow_rate,
                        double burn_start_time, double burn_duration, 
                        double burn_heading, double burn_throttle);


void createCSV(char* FILENAME);
void readCSV(char* FILENAME, body_properties_t** gb, int* num_bodies);
void readSpacecraftCSV(char* FILENAME, spacecraft_properties_t** sc, int* num_craft);
void readBodyJSON(char* FILENAME, body_properties_t** gb, int* num_bodies);
void readSpacecraftJSON(char* FILENAME, spacecraft_properties_t** sc, int* num_craft);

void runCalculations(body_properties_t** gb, spacecraft_properties_t** sc, window_params_t* wp, int num_bodies, int num_craft); // MAIN CALCULATION LOOP -- DOES ALL THE MATH

#endif