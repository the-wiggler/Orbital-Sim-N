#ifndef SIMULATION_H
#define SIMULATION_H

#include "../types.h"

double calculateTotalSystemEnergy(const sim_properties_t* sim);
void resetSim(sim_properties_t* sim);
void updateSystemOrbitalElements(sim_properties_t* sim);
void runCalculations(sim_properties_t* sim);
void calculateOrbitalElements(orbital_elements_t* target_orbital_elements, const vec3* target_pos, const vec3* target_vel, const body_t* orbiting_body);
void cleanup(const sim_properties_t* sim);

#endif
