#ifndef BODIES_H
#define BODIES_H

#include "../types.h"

void body_calculateGravForce(sim_properties_t* sim, int i, int j);
void body_updateMotion(body_t* body, double dt);
void body_updateRotation(body_t* body, double dt);
void body_calculateKineticEnergy(body_t* body);
void body_addOrbitalBody(body_properties_t* gb, const char* name, double mass, double radius, vec3 pos, vec3 vel);

#endif
