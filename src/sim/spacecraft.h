#ifndef SPACECRAFT_H
#define SPACECRAFT_H

#include "../types.h"

typedef struct {
    vec3 optimal_delta_v;
    int steps_until_burn;
    double orbit_time_step;
    double expected_time_of_contact;
} delta_v_optimizer_data;

vec3 craft_calculateGravForce(const sim_properties_t* sim, int craft_idx, int body_idx);
vec3 craft_calculateJ2Force(const sim_properties_t* sim, int craft_idx, int body_idx);
void craft_updateMotion(spacecraft_t* craft, double delta_t);
void craft_applyThrust(spacecraft_t* craft);
void craft_checkBurnSchedule(spacecraft_t* craft, const body_properties_t* global_bodies, double sim_time);
void craft_consumeFuel(spacecraft_t* craft, double delta_t);
delta_v_optimizer_data craft_findOptimalDeltaV(const sim_properties_t* sim, int craft_id, int orbit_sample_points, int point_samples);
burn_properties_t craft_autoDeltaVOptimization(const sim_properties_t* sim, int craft_id, int orbit_sample_points, int point_samples);
void craft_addSpacecraft(spacecraft_properties_t* global_spacecraft, const char* name,
                        vec3 pos, vec3 vel,
                        double dry_mass, double fuel_mass, double thrust,
                        double specific_impulse, double mass_flow_rate,
                        double attitude, double moment_of_inertia,
                        double nozzle_gimbal_range, double nozzle_velocity,
                        const burn_properties_t* burns, int num_burns);
void craft_findClosestPlanet(spacecraft_t* craft, const body_properties_t* global_bodies);

#endif
