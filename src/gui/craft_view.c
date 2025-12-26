//
// Created by java on 12/11/2025.
//
#include "../globals.h"
#include "../types.h"
#include "craft_view.h"
#include <math.h>
#include <stdio.h>

#include "SDL_engine.h"

// calculates the distance and heading of a craft to another body
typedef struct {
    double heading;
    double distance_x, distance_y, distance;
    double rel_v_x, rel_v_y, rel_v;
} rel_body_vector_t;

rel_body_vector_t craft_calcVectorToBody(const spacecraft_properties_t* sc, const body_properties_t* gb, const int craft_id, const int body_id) {
    // relative distance
    const double dist_x = gb->pos_x[body_id] - sc->pos_x[craft_id];
    const double dist_y = gb->pos_y[body_id] - sc->pos_y[craft_id];
    const double dist = sqrt(pow(dist_x, 2) + pow(dist_y, 2));

    // relative velocity
    const double rvx = gb->vel_x[body_id] - sc->vel_x[craft_id];
    const double rvy = gb->vel_y[body_id] - sc->vel_y[craft_id];
    const double rv = sqrt(pow(rvx, 2) + pow(rvy, 2));

    // heading (radians)
    const double hd = atan2(-dist_y, dist_x);

    rel_body_vector_t result;
    result.heading = hd;
    result.distance_x = dist_x;
    result.distance_y = dist_y;
    result.distance = dist;
    result.rel_v_x = rvx;
    result.rel_v_y = rvy;
    result.rel_v = rv;
    return result;
}

void craft_RenderCraftView(sim_properties_t* sim) {
    (void)sim;  // Unused parameter
    // TODO: Implement OpenGL rendering for craft view
}
