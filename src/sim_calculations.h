#ifndef SIM_CALCULATIONS_H
#define SIM_CALCULATIONS_H

#include <SDL3_ttf/SDL_ttf.h>
#include "config.h"

extern const double G;
extern TTF_Font* g_font;

// orbital body stuff
void body_calculateGravForce(const body_properties_t* bodies, int i, int j, window_params_t* wp);
void body_updateMotion(const body_properties_t* bodies, int i, double dt);
void body_transformCoordinates(const body_properties_t* bodies, int i, window_params_t window_params);
void body_calculateKineticEnergy(const body_properties_t* bodies, int i);
float body_calculateVisualRadius(const body_properties_t* bodies, int i, window_params_t wp);
void body_addOrbitalBody(body_properties_t* gb, const char* name, double mass, double radius, double x_pos, double y_pos, double x_vel, double y_vel);
void resetSim(window_params_t* wp, body_properties_t* gb, spacecraft_properties_t* sc);

// system energy calculation
double calculateTotalSystemEnergy(const body_properties_t* gb, const spacecraft_properties_t* sc);

// spacecraft body stuff
void craft_calculateGravForce(const spacecraft_properties_t* sc, int i, const body_properties_t* bodies, int j, window_params_t* wp);
void craft_transformCoordinates(const spacecraft_properties_t* sc, int i, window_params_t wp);
void craft_updateMotion(const spacecraft_properties_t* sc, int i, double dt);
void craft_applyThrust(const spacecraft_properties_t* sc, int i);
void craft_consumeFuel(const spacecraft_properties_t* sc, int i, double dt);
void craft_addSpacecraft(spacecraft_properties_t* sc, const char* name,
                        double x_pos, double y_pos, double x_vel, double y_vel,
                        double dry_mass, double fuel_mass, double thrust,
                        double specific_impulse, double mass_flow_rate,
                        double attitude, double moment_of_inertia,
                        double nozzle_gimbal_range,
                        double burn_start_time, double burn_duration,
                        double burn_heading, double burn_throttle);

void runCalculations(const body_properties_t* gb, const spacecraft_properties_t* sc, window_params_t* wp); // MAIN CALCULATION LOOP -- DOES ALL THE MATH

#endif