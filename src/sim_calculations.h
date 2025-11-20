#ifndef SIM_CALCULATIONS_H
#define SIM_CALCULATIONS_H

#include <math.h>
#include <SDL3/SDL.h>
#include <SDL3_ttf/SDL_ttf.h>
#include <stdbool.h>
#include "config.h"

extern const double G;
extern TTF_Font* g_font;

void calculateForce(body_properties_t* b, body_properties_t b2);
void updateMotion(body_properties_t* b, double dt);
void transformCoordinates(body_properties_t* b, window_params_t window_params);
void calculateKineticEnergy(body_properties_t* b);
double calculatePotentialEnergy(body_properties_t b, const char* target_name, body_properties_t* gb, int num_bodies);
int calculateVisualRadius(body_properties_t* body, window_params_t wp);
void addOrbitalBody(body_properties_t** gb, int* num_bodies, char* name, double mass, double x_pos, double y_pos, double x_vel, double y_vel);
void resetSim(double* sim_time, body_properties_t** gb, int* num_bodies);

void createCSV(char* FILENAME);
void readCSV(char* FILENAME, body_properties_t** gb, int* num_bodies);
void exportOrbitalEnergyData(char* FILENAME, body_properties_t* gb, int num_bodies, double sim_time);

void runCalculations(body_properties_t** gb, window_params_t* wp, int num_bodies); // MAIN CALCULATION LOOP -- DOES ALL THE MATH

#endif