#ifndef BODIES_H
#define BODIES_H

#include "../types.h"

vec3 body_calculateGravForce(const sim_properties_t* sim, int force_recipient_id, int force_applier_id);
// TODO: edit these function so that they return types instead of void
void body_updateRotation(body_t* body, double delta_t);
void body_precomputeRotations(body_properties_t* global_bodies, double delta_t);
void body_calculateSOI(body_properties_t* global_bodies);
void body_findClosestPlanet(body_t* body, const body_properties_t* global_bodies);
void body_addOrbitalBody(body_properties_t* global_bodies, const char* name, double mass, double radius, vec3 pos, vec3 vel);

#endif
