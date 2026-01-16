#ifndef SPACECRAFT_H
#define SPACECRAFT_H

#include "../types.h"

void craft_calculateGravForce(sim_properties_t* sim, int craft_idx, int body_idx);
void craft_calculateOrbitalElements(spacecraft_t* craft, const body_t* body);
void craft_updateMotion(spacecraft_t* craft, double dt);
void craft_applyThrust(spacecraft_t* craft);
void craft_checkBurnSchedule(spacecraft_t* craft, const body_properties_t* gb, double sim_time);
void craft_consumeFuel(spacecraft_t* craft, double dt);
void craft_addSpacecraft(spacecraft_properties_t* gs, const char* name,
                        vec3 pos, vec3 vel,
                        double dry_mass, double fuel_mass, double thrust,
                        double specific_impulse, double mass_flow_rate,
                        double attitude, double moment_of_inertia,
                        double nozzle_gimbal_range,
                        const burn_properties_t* burns, int num_burns);
void craft_findClosestPlanet(spacecraft_t* craft, body_properties_t* gb);
#endif
