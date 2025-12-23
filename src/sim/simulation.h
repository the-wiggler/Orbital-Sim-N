#ifndef SIMULATION_H
#define SIMULATION_H

#include "../types.h"

double calculateTotalSystemEnergy(const body_properties_t* gb, const spacecraft_properties_t* sc);
void resetSim(window_params_t* wp, body_properties_t* gb, spacecraft_properties_t* sc);
void runCalculations(const body_properties_t* gb, const spacecraft_properties_t* sc, window_params_t* wp);
void cleanup(void* args);

#endif
