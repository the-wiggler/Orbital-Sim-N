#ifndef SPACECRAFT_H
#define SPACECRAFT_H

#include "../types.h"

void craft_calculateGravForce(sim_properties_t* sim, int craft_idx, int body_idx);
void craft_transformCoordinates(const spacecraft_properties_t* sc, int i, window_params_t wp);
void craft_updateMotion(const spacecraft_properties_t* sc, int i, double dt);
void craft_applyThrust(const spacecraft_properties_t* sc, int i);
void craft_checkBurnSchedule(const spacecraft_properties_t* sc, int i, const body_properties_t* gb, double sim_time);
void craft_consumeFuel(const spacecraft_properties_t* sc, int i, double dt);
void craft_addSpacecraft(spacecraft_properties_t* sc, const char* name,
                        double x_pos, double y_pos, double z_pos, double x_vel, double y_vel, double z_vel,
                        double dry_mass, double fuel_mass, double thrust,
                        double specific_impulse, double mass_flow_rate,
                        double attitude, double moment_of_inertia,
                        double nozzle_gimbal_range,
                        const burn_properties_t* burns, int num_burns);

#endif
