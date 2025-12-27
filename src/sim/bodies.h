#ifndef BODIES_H
#define BODIES_H

#include "../types.h"

void body_calculateGravForce(sim_properties_t* sim, int i, int j);
void body_updateMotion(const body_properties_t* bodies, int i, double dt);
void body_transformCoordinates(const body_properties_t* bodies, int i, window_params_t window_params);
void body_calculateKineticEnergy(const body_properties_t* bodies, int i);
float body_calculateVisualRadius(const body_properties_t* bodies, int i, window_params_t wp);
void body_addOrbitalBody(body_properties_t* gb, const char* name, double mass, double radius, double x_pos, double y_pos, double z_pos, double x_vel, double y_vel, double z_vel);

#endif
