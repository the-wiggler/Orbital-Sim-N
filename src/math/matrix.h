//
// Created by java on 12/27/25.
//

#ifndef ORBITSIMULATION_MATRIX_H
#define ORBITSIMULATION_MATRIX_H

#include "../types.h"
#include <math.h>

// normalize a 3d vector
static inline void normalize_3d(float v[3]) {
    float length = sqrtf(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    if (length > 0.0f) {
        v[0] /= length;
        v[1] /= length;
        v[2] /= length;
    }
}

// cross product of two 3d vectors (float)
static inline void cross_product_3d(const float a[3], const float b[3], float result[3]) {
    result[0] = a[1] * b[2] - a[2] * b[1];
    result[1] = a[2] * b[0] - a[0] * b[2];
    result[2] = a[0] * b[1] - a[1] * b[0];
}

// cross product of two vec3 vectors (double precision)
static inline vec3 cross_product_vec3(const vec3 a, const vec3 b) {
    vec3 result;
    result.x = a.y * b.z - a.z * b.y;
    result.y = a.z * b.x - a.x * b.z;
    result.z = a.x * b.y - a.y * b.x;
    return result;
}

static inline vec3 vec3_zero(void) {
    return (vec3){0.0, 0.0, 0.0};
}

// add
static inline vec3 vec3_add(vec3 a, vec3 b) {
    return (vec3){a.x + b.x, a.y + b.y, a.z + b.z};
}

// subtract
static inline vec3 vec3_sub(vec3 a, vec3 b) {
    return (vec3){a.x - b.x, a.y - b.y, a.z - b.z};
}

// change magnitude of a vector
static inline vec3 vec3_scale(vec3 v, double s) {
    return (vec3){v.x * s, v.y * s, v.z * s};
}

// dot product
static inline double vec3_dot(vec3 a, vec3 b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

// cross product
static inline vec3 vec3_cross(vec3 a, vec3 b) {
    return (vec3){a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}

// vector division by scalar
static inline vec3 vec3_scalar_div(vec3 a, double s) {
    return (vec3){a.x / s, a.y / s, a.z / s};
}

// squared magnitude
static inline double vec3_mag_sq(vec3 v) {
    return v.x * v.x + v.y * v.y + v.z * v.z;
}

// magnitude of vector
static inline double vec3_mag(vec3 v) {
    return sqrt(vec3_mag_sq(v));
}

// normalize
static inline vec3 vec3_normalize(vec3 v) {
    double mag = vec3_mag(v);
    if (mag > 0.0) {
        return vec3_scale(v, 1.0 / mag);
    }
    return v;
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
static inline vec3_f mat4_transformPoint(mat4 m, vec3_f point) {
    vec3_f result;
    result.x = m.m[0] * point.x + m.m[4] * point.y + m.m[8]  * point.z + m.m[12];
    result.y = m.m[1] * point.x + m.m[5] * point.y + m.m[9]  * point.z + m.m[13];
    result.z = m.m[2] * point.x + m.m[6] * point.y + m.m[10] * point.z + m.m[14];
    return result;
}

// quaternion math
static inline quaternion_t quaternionMul(quaternion_t q1, quaternion_t q2) {
    quaternion_t r;
    r.w = q1.w * q2.w - q1.x * q2.x - q1.y * q2.y - q1.z * q2.z;
    r.x = q1.w * q2.x + q1.x * q2.w + q1.y * q2.z - q1.z * q2.y;
    r.y = q1.w * q2.y - q1.x * q2.z + q1.y * q2.w + q1.z * q2.x;
    r.z = q1.w * q2.z + q1.x * q2.y - q1.y * q2.x + q1.z * q2.w;
    return r;
}

static inline quaternion_t quaternionFromAxisAngle(const vec3 axis, double angle) {
    quaternion_t q;
    double half_angle = angle / 2.0;
    double s = sin(half_angle);
    double norm = sqrt(axis.x*axis.x + axis.y*axis.y + axis.z*axis.z);
    q.w = cos(half_angle);
    q.x = axis.x / norm * s;
    q.y = axis.y / norm * s;
    q.z = axis.z / norm * s;
    return q;
}

static inline vec3 quaternionRotate(const quaternion_t q, const vec3 v) {
    vec3 r;
    quaternion_t qv = {0.0, v.x, v.y, v.z};
    quaternion_t q_conj = {q.w, -q.x, -q.y, -q.z};
    quaternion_t temp = quaternionMul(q, qv);
    quaternion_t rq = quaternionMul(temp, q_conj);
    r.x = rq.x;
    r.y = rq.y;
    r.z = rq.z;
    return r;
}

// creates a quaternion that rotates from vector 'from' to vector 'to'
// both vectors should be normalized
static inline quaternion_t quaternionFromTwoVectors(vec3 from, vec3 to) {
    // normalize the input vectors
    double from_len = sqrt(from.x * from.x + from.y * from.y + from.z * from.z);
    double to_len = sqrt(to.x * to.x + to.y * to.y + to.z * to.z);

    if (from_len > 0.0) {
        from.x /= from_len;
        from.y /= from_len;
        from.z /= from_len;
    }

    if (to_len > 0.0) {
        to.x /= to_len;
        to.y /= to_len;
        to.z /= to_len;
    }

    // compute dot product (cosine of angle)
    double dot = from.x * to.x + from.y * to.y + from.z * to.z;

    // check if vectors are already aligned
    if (dot > 0.999999) {
        // vectors are parallel - return identity quaternion
        return (quaternion_t){1.0, 0.0, 0.0, 0.0};
    }

    // check if vectors are opposite
    if (dot < -0.999999) {
        // vectors are anti-parallel - rotate 180 degrees around any perpendicular axis
        // find a perpendicular axis
        vec3 axis;
        if (fabs(from.x) < 0.9) {
            axis = (vec3){1.0, 0.0, 0.0};
        } else {
            axis = (vec3){0.0, 1.0, 0.0};
        }

        // cross product to get perpendicular
        vec3 perp;
        perp.x = from.y * axis.z - from.z * axis.y;
        perp.y = from.z * axis.x - from.x * axis.z;
        perp.z = from.x * axis.y - from.y * axis.x;

        // normalize
        double perp_len = sqrt(perp.x * perp.x + perp.y * perp.y + perp.z * perp.z);
        perp.x /= perp_len;
        perp.y /= perp_len;
        perp.z /= perp_len;

        // 180 degree rotation around perpendicular axis
        return (quaternion_t){0.0, perp.x, perp.y, perp.z};
    }

    // compute rotation axis (cross product)
    vec3 axis;
    axis.x = from.y * to.z - from.z * to.y;
    axis.y = from.z * to.x - from.x * to.z;
    axis.z = from.x * to.y - from.y * to.x;

    // quaternion components
    quaternion_t q;
    q.w = sqrt((1.0 + dot) * 0.5);
    double s = sqrt((1.0 - dot) * 0.5);
    q.x = axis.x * s;
    q.y = axis.y * s;
    q.z = axis.z * s;

    return q;
}

static inline mat4 quaternionToMatrix(const quaternion_t q) {
    float xx = (float)(q.x * q.x);
    float yy = (float)(q.y * q.y);
    float zz = (float)(q.z * q.z);
    float xy = (float)(q.x * q.y);
    float xz = (float)(q.x * q.z);
    float yz = (float)(q.y * q.z);
    float wx = (float)(q.w * q.x);
    float wy = (float)(q.w * q.y);
    float wz = (float)(q.w * q.z);

    mat4 mat;
    // column 0
    mat.m[0] = 1.0f - 2.0f * (yy + zz);
    mat.m[1] = 2.0f * (xy + wz);
    mat.m[2] = 2.0f * (xz - wy);
    mat.m[3] = 0.0f;

    // column 1
    mat.m[4] = 2.0f * (xy - wz);
    mat.m[5] = 1.0f - 2.0f * (xx + zz);
    mat.m[6] = 2.0f * (yz + wx);
    mat.m[7] = 0.0f;

    // column 2
    mat.m[8] = 2.0f * (xz + wy);
    mat.m[9] = 2.0f * (yz - wx);
    mat.m[10] = 1.0f - 2.0f * (xx + yy);
    mat.m[11] = 0.0f;

    // column 3
    mat.m[12] = 0.0f;
    mat.m[13] = 0.0f;
    mat.m[14] = 0.0f;
    mat.m[15] = 1.0f;

    return mat;
}

#endif //ORBITSIMULATION_MATRIX_H