//
// Created by java on 12/27/25.
//

#ifndef ORBITSIMULATION_MATRIX_H
#define ORBITSIMULATION_MATRIX_H

#include "../types.h"
#include <math.h>

// normalize a 3d vector
static inline void normalize_3d(float vec[3]) {
    const float length = sqrtf((vec[0] * vec[0]) + (vec[1] * vec[1]) + (vec[2] * vec[2]));
    if (length > 0.0F) {
        vec[0] /= length;
        vec[1] /= length;
        vec[2] /= length;
    }
}

// cross product of two 3d vectors (float)
static inline void cross_product_3d(const float arg1[3], const float arg2[3], float result[3]) {
    result[0] = (arg1[1] * arg2[2]) - (arg1[2] * arg2[1]);
    result[1] = (arg1[2] * arg2[0]) - (arg1[0] * arg2[2]);
    result[2] = (arg1[0] * arg2[1]) - (arg1[1] * arg2[0]);
}

// cross product of two vec3 vectors (double precision)
static inline vec3 cross_product_vec3(const vec3 vec_a, const vec3 vec_b) {
    vec3 result;
    result.x = (vec_a.y * vec_b.z) - (vec_a.z * vec_b.y);
    result.y = (vec_a.z * vec_b.x) - (vec_a.x * vec_b.z);
    result.z = (vec_a.x * vec_b.y) - (vec_a.y * vec_b.x);
    return result;
}

static inline vec3 vec3_zero(void) {
    return (vec3){0.0, 0.0, 0.0};
}

// add
static inline vec3 vec3_add(const vec3 vec_a, const vec3 vec_b) {
    return (vec3){vec_a.x + vec_b.x, vec_a.y + vec_b.y, vec_a.z + vec_b.z};
}

// subtract
static inline vec3 vec3_sub(const vec3 vec_a, const vec3 vec_b) {
    return (vec3){vec_a.x - vec_b.x, vec_a.y - vec_b.y, vec_a.z - vec_b.z};
}
static inline vec3_f vec3_f_sub(const vec3_f vec_a, const vec3_f vec_b) {
    return (vec3_f){vec_a.x - vec_b.x, vec_a.y - vec_b.y, vec_a.z - vec_b.z};
}

// change magnitude of a vector
static inline vec3 vec3_scale(const vec3 vec, const double scalar) {
    return (vec3){vec.x * scalar, vec.y * scalar, vec.z * scalar};
}

// dot product
static inline double vec3_dot(const vec3 vec_a, const vec3 vec_b) {
    return (vec_a.x * vec_b.x) + (vec_a.y * vec_b.y) + (vec_a.z * vec_b.z);
}

// cross product
static inline vec3 vec3_cross(const vec3 vec_a, const vec3 vec_b) {
    return (vec3){(vec_a.y * vec_b.z) - (vec_a.z * vec_b.y), (vec_a.z * vec_b.x) - (vec_a.x * vec_b.z), (vec_a.x * vec_b.y) - (vec_a.y * vec_b.x)};
}

// vector division by scalar
static inline vec3 vec3_scalar_div(const vec3 vec, const double scalar) {
    return (vec3){vec.x / scalar, vec.y / scalar, vec.z / scalar};
}

// squared magnitude
static inline double vec3_mag_sq(const vec3 vec) {
    return (vec.x * vec.x) + (vec.y * vec.y) + (vec.z * vec.z);
}
static inline float vec3_f_mag_sq(const vec3_f vec) {
    return (vec.x * vec.x) + (vec.y * vec.y) + (vec.z * vec.z);
}

// magnitude of vector
static inline double vec3_mag(const vec3 vec) {
    return sqrt(vec3_mag_sq(vec));
}
static inline float vec3_f_mag(const vec3_f vec) {
    return sqrtf(vec3_f_mag_sq(vec));
}

// normalize
static inline vec3 vec3_normalize(const vec3 vec) {
    const double mag = vec3_mag(vec);
    if (mag > 0.0) {
        return vec3_scale(vec, 1.0 / mag);
    }
    return vec;
}
static inline vec3_f vec3_f_normalize(const vec3_f vec) {
    const float mag = vec3_f_mag(vec);
    if (mag > 0.0F) {
        return (vec3_f){vec.x / mag, vec.y / mag, vec.z / mag};
    }
    return vec;
}

// cross product for vec3_f
static inline vec3_f vec3_f_cross(const vec3_f vec_a, const vec3_f vec_b) {
    return (vec3_f){
        (vec_a.y * vec_b.z) - (vec_a.z * vec_b.y),
        (vec_a.z * vec_b.x) - (vec_a.x * vec_b.z),
        (vec_a.x * vec_b.y) - (vec_a.y * vec_b.x)
    };
}

// dot product for vec3_f
static inline float vec3_f_dot(const vec3_f vec_a, const vec3_f vec_b) {
    return (vec_a.x * vec_b.x) + (vec_a.y * vec_b.y) + (vec_a.z * vec_b.z);
}

// creates an identity matrix
static inline mat4 mat4_identity(void) {
    const mat4 matrix = {.m = {
        1, 0, 0, 0,
        0, 1, 0, 0,
        0, 0, 1, 0,
        0, 0, 0, 1
    }};
    return matrix;
}

// creates a translation matrix
static inline mat4 mat4_translation(const float pos_x, const float pos_y, const float pos_z) {
    const mat4 matrix = {.m = {
        1, 0, 0, 0,
        0, 1, 0, 0,
        0, 0, 1, 0,
        pos_x, pos_y, pos_z, 1
    }};
    return matrix;
}

// creates a scale matrix
static inline mat4 mat4_scale(const float scale_x, const float scale_y, const float scale_z) {
    const mat4 matrix = {.m = {
        scale_x, 0,       0,       0,
        0,       scale_y, 0,       0,
        0,       0,       scale_z, 0,
        0,       0,       0,       1
    }};
    return matrix;
}

// creates a rotation matrix around the X axis
static inline mat4 mat4_rotationX(const float angle) {
    const float cosine = cosf(angle);
    const float sine = sinf(angle);
    const mat4 matrix = {.m = {
        1, 0,       0,        0,
        0, cosine, -sine,     0,
        0, sine,    cosine,   0,
        0, 0,       0,        1
    }};
    return matrix;
}

// creates a rotation matrix around the Y axis
static inline mat4 mat4_rotationY(const float angle) {
    const float cosine = cosf(angle);
    const float sine = sinf(angle);
    const mat4 matrix = {.m = {
         cosine, 0, sine,    0,
         0,      1, 0,       0,
        -sine,   0, cosine,  0,
         0,      0, 0,       1
    }};
    return matrix;
}

// creates a rotation matrix around the Z axis
static inline mat4 mat4_rotationZ(const float angle) {
    const float cosine = cosf(angle);
    const float sine = sinf(angle);
    const mat4 matrix = {.m = {
        cosine, -sine,   0, 0,
        sine,    cosine, 0, 0,
        0,       0,      1, 0,
        0,       0,      0, 1
    }};
    return matrix;
}

// matrix multiplication
static inline mat4 mat4_mul(const mat4 arg1, const mat4 arg2) {
    mat4 result;

    for (int col = 0; col < 4; ++col) {
        for (int row = 0; row < 4; ++row) {
            result.m[(col*4) + row] =
                (arg1.m[(0*4) + row] * arg2.m[(col*4) + 0]) +
                (arg1.m[(1*4) + row] * arg2.m[(col*4) + 1]) +
                (arg1.m[(2*4) + row] * arg2.m[(col*4) + 2]) +
                (arg1.m[(3*4) + row] * arg2.m[(col*4) + 3]);
        }
    }
    return result;
}

// transform a point by a matrix (for rotating camera position)
static inline vec3_f vec3_transformByMat4(const mat4 matrix, const vec3_f point) {
    vec3_f result;
    result.x = (matrix.m[0] * point.x) + (matrix.m[4] * point.y) + (matrix.m[8]  * point.z) + matrix.m[12];
    result.y = (matrix.m[1] * point.x) + (matrix.m[5] * point.y) + (matrix.m[9]  * point.z) + matrix.m[13];
    result.z = (matrix.m[2] * point.x) + (matrix.m[6] * point.y) + (matrix.m[10] * point.z) + matrix.m[14];
    return result;
}

// quaternion math
static inline quaternion_t quaternionMul(const quaternion_t quat1, const quaternion_t quat2) {
    quaternion_t result;
    result.w = (quat1.w * quat2.w) - (quat1.x * quat2.x) - (quat1.y * quat2.y) - (quat1.z * quat2.z);
    result.x = (quat1.w * quat2.x) + (quat1.x * quat2.w) + (quat1.y * quat2.z) - (quat1.z * quat2.y);
    result.y = (quat1.w * quat2.y) - (quat1.x * quat2.z) + (quat1.y * quat2.w) + (quat1.z * quat2.x);
    result.z = (quat1.w * quat2.z) + (quat1.x * quat2.y) - (quat1.y * quat2.x) + (quat1.z * quat2.w);
    return result;
}

static inline quaternion_t quaternionFromAxisAngle(const vec3 axis, const double angle) {
    quaternion_t quaternion;
    const double half_angle = angle / 2.0;
    const double sine = sin(half_angle);
    const double norm = sqrt((axis.x*axis.x) + (axis.y*axis.y) + (axis.z*axis.z));
    quaternion.w = cos(half_angle);
    quaternion.x = (axis.x / norm) * sine;
    quaternion.y = (axis.y / norm) * sine;
    quaternion.z = (axis.z / norm) * sine;
    return quaternion;
}

static inline vec3 quaternionRotate(const quaternion_t quaternion, const vec3 vec) {
    vec3 result;
    const quaternion_t quat_vec = {0.0, vec.x, vec.y, vec.z};
    const quaternion_t q_conj = {quaternion.w, -quaternion.x, -quaternion.y, -quaternion.z};
    const quaternion_t temp = quaternionMul(quaternion, quat_vec);
    const quaternion_t result_quat = quaternionMul(temp, q_conj);
    result.x = result_quat.x;
    result.y = result_quat.y;
    result.z = result_quat.z;
    return result;
}

// creates a quaternion that rotates from vector 'from' to vector 'target_vec'
// both vectors should be normalized
static inline quaternion_t quaternionFromTwoVectors(vec3 from, vec3 target_vec) {
    // normalize the input vectors
    const double from_len = sqrt((from.x * from.x) + (from.y * from.y) + (from.z * from.z));
    const double target_len = sqrt((target_vec.x * target_vec.x) + (target_vec.y * target_vec.y) + (target_vec.z * target_vec.z));

    if (from_len > 0.0) {
        from.x /= from_len;
        from.y /= from_len;
        from.z /= from_len;
    }

    if (target_len > 0.0) {
        target_vec.x /= target_len;
        target_vec.y /= target_len;
        target_vec.z /= target_len;
    }

    // compute dot product (cosine of angle)
    const double dot = (from.x * target_vec.x) + (from.y * target_vec.y) + (from.z * target_vec.z);

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
        const double perp_len = sqrt((perp.x * perp.x) + (perp.y * perp.y) + (perp.z * perp.z));
        perp.x /= perp_len;
        perp.y /= perp_len;
        perp.z /= perp_len;

        // 180 degree rotation around perpendicular axis
        return (quaternion_t){0.0, perp.x, perp.y, perp.z};
    }

    // compute rotation axis (cross product)
    vec3 axis;
    axis.x = from.y * target_vec.z - from.z * target_vec.y;
    axis.y = from.z * target_vec.x - from.x * target_vec.z;
    axis.z = from.x * target_vec.y - from.y * target_vec.x;

    // normalize the axis (cross product magnitude is sin(angle), not 1)
    const double axis_len = sqrt((axis.x * axis.x) + (axis.y * axis.y) + (axis.z * axis.z));
    if (axis_len > 0.0) {
        axis.x /= axis_len;
        axis.y /= axis_len;
        axis.z /= axis_len;
    }

    // quaternion components
    quaternion_t quat;
    quat.w = sqrt((1.0 + dot) * 0.5);
    const double sine = sqrt((1.0 - dot) * 0.5);
    quat.x = axis.x * sine;
    quat.y = axis.y * sine;
    quat.z = axis.z * sine;

    return quat;
}

static inline mat4 quaternionToMatrix(const quaternion_t quaternion) {
    const float quat_xx = (float)(quaternion.x * quaternion.x);
    const float quat_yy = (float)(quaternion.y * quaternion.y);
    const float quat_zz = (float)(quaternion.z * quaternion.z);
    const float quat_xy = (float)(quaternion.x * quaternion.y);
    const float quat_xz = (float)(quaternion.x * quaternion.z);
    const float quat_yz = (float)(quaternion.y * quaternion.z);
    const float quat_wx = (float)(quaternion.w * quaternion.x);
    const float quat_wy = (float)(quaternion.w * quaternion.y);
    const float quat_wz = (float)(quaternion.w * quaternion.z);

    mat4 matrix;
    // column 0
    matrix.m[0] = 1.0F - (2.0F * (quat_yy + quat_zz));
    matrix.m[1] = 2.0F * (quat_xy + quat_wz);
    matrix.m[2] = 2.0F * (quat_xz - quat_wy);
    matrix.m[3] = 0.0F;

    // column 1
    matrix.m[4] = 2.0F * (quat_xy - quat_wz);
    matrix.m[5] = 1.0F - (2.0F * (quat_xx + quat_zz));
    matrix.m[6] = 2.0F * (quat_yz + quat_wx);
    matrix.m[7] = 0.0F;

    // column 2
    matrix.m[8] = 2.0F * (quat_xz + quat_wy);
    matrix.m[9] = 2.0F * (quat_yz - quat_wx);
    matrix.m[10] = 1.0F - (2.0F * (quat_xx + quat_yy));
    matrix.m[11] = 0.0F;

    // column 3
    matrix.m[12] = 0.0F;
    matrix.m[13] = 0.0F;
    matrix.m[14] = 0.0F;
    matrix.m[15] = 1.0F;

    return matrix;
}

#endif //ORBITSIMULATION_MATRIX_H