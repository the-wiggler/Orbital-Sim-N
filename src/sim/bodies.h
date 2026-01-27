#ifndef BODIES_H
#define BODIES_H

#include "../types.h"

void body_calculateGravForce(sim_properties_t* sim, int force_recipient, int force_applier);
void body_updateRotation(body_t* body, double dt);
void body_calculateKineticEnergy(body_t* body);
void body_calculateSOI(body_properties_t* gb);
void body_addOrbitalBody(body_properties_t* gb, const char* name, double mass, double radius, vec3 pos, vec3 vel);

#endif
