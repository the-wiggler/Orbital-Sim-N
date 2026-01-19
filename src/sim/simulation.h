#ifndef SIMULATION_H
#define SIMULATION_H

#include "../types.h"

double calculateTotalSystemEnergy(const sim_properties_t* sim);
void resetSim(sim_properties_t* sim);
void runCalculations(sim_properties_t* sim);
void cleanup(const sim_properties_t* sim);

#endif
