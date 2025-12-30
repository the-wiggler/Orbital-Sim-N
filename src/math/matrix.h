//
// Created by java on 12/27/25.
//

#ifndef ORBITSIMULATION_MATRIX_H
#define ORBITSIMULATION_MATRIX_H

#include "../types.h"
#include <math.h>

// normalize a 3D vector
static inline void normalize_3d(float v[3]) {
    float length = sqrtf(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    if (length > 0.0f) {
        v[0] /= length;
        v[1] /= length;
        v[2] /= length;
    }
}

// cross product of two 3D vectors
static inline void cross_product_3d(const float a[3], const float b[3], float result[3]) {
    result[0] = a[1] * b[2] - a[2] * b[1];
    result[1] = a[2] * b[0] - a[0] * b[2];
    result[2] = a[0] * b[1] - a[1] * b[0];
}

// creates an identity matrix
static inline mat4 mat4_identity(void) {
    mat4 m = {.m = {
        1, 0, 0, 0,
        0, 1, 0, 0,
        0, 0, 1, 0,
        0, 0, 0, 1
    }};
    return m;
}

// creates a translation matrix
static inline mat4 mat4_translation(float x, float y, float z) {
    mat4 m = {.m = {
        1, 0, 0, 0,
        0, 1, 0, 0,
        0, 0, 1, 0,
        x, y, z, 1
    }};
    return m;
}

// creates a scale matrix
static inline mat4 mat4_scale(float sx, float sy, float sz) {
    mat4 m = {.m = {
        sx, 0,  0,  0,
        0,  sy, 0,  0,
        0,  0,  sz, 0,
        0,  0,  0,  1
    }};
    return m;
}

// creates a rotation matrix around the X axis
static inline mat4 mat4_rotationX(float angle) {
    float c = cosf(angle);
    float s = sinf(angle);
    mat4 m = {.m = {
        1, 0,  0, 0,
        0, c, -s, 0,
        0, s,  c, 0,
        0, 0,  0, 1
    }};
    return m;
}

// creates a rotation matrix around the Y axis
static inline mat4 mat4_rotationY(float angle) {
    float c = cosf(angle);
    float s = sinf(angle);
    mat4 m = {.m = {
         c, 0, s, 0,
         0, 1, 0, 0,
        -s, 0, c, 0,
         0, 0, 0, 1
    }};
    return m;
}

// creates a rotation matrix around the Z axis
static inline mat4 mat4_rotationZ(float angle) {
    float c = cosf(angle);
    float s = sinf(angle);
    mat4 m = {.m = {
        c, -s, 0, 0,
        s,  c, 0, 0,
        0,  0, 1, 0,
        0,  0, 0, 1
    }};
    return m;
}

// matrix multiplication
static inline mat4 mat4_mul(mat4 a, mat4 b) {
    mat4 r;

    for (int col = 0; col < 4; ++col) {
        for (int row = 0; row < 4; ++row) {
            r.m[col*4 + row] =
                a.m[0*4 + row] * b.m[col*4 + 0] +
                a.m[1*4 + row] * b.m[col*4 + 1] +
                a.m[2*4 + row] * b.m[col*4 + 2] +
                a.m[3*4 + row] * b.m[col*4 + 3];
        }
    }
    return r;
}

// transform a point by a matrix (for rotating camera position)
static inline coord_t mat4_transformPoint(mat4 m, coord_t point) {
    coord_t result;
    result.x = m.m[0] * point.x + m.m[4] * point.y + m.m[8]  * point.z + m.m[12];
    result.y = m.m[1] * point.x + m.m[5] * point.y + m.m[9]  * point.z + m.m[13];
    result.z = m.m[2] * point.x + m.m[6] * point.y + m.m[10] * point.z + m.m[14];
    return result;
}

#endif //ORBITSIMULATION_MATRIX_H