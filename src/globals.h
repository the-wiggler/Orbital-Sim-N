#ifndef GLOBALS_H
#define GLOBALS_H

//////////////////////
// COMPILE CONFIG:  //
//////////////////////

static const double G_C = 6.67430E-11;
static const double NEAR_ZERO = 1e-10;  // threshold for detecting near-zero values in numerical calculations
static const double STANDARD_GRAVITY = 9.80665;

static const double N_PI = 3.14159265358979323846;
static const float PI_f = 3.14159265358979323846F;

static const double PI_OVER_6    =  0.52359877559829877312;  // π/6
static const double PI_OVER_4    =  0.78539816339744830962;  // π/4
static const double PI_OVER_3    =  1.04719755119659763133;  // π/3
static const double PI_OVER_2    =  1.57079632679489661923;  // π/2
static const double TWO_PI       =  6.28318530717958623200;  // 2π

static const float PI_OVER_6_f   =  0.5235987756F;  // π/6
static const float PI_OVER_4_f   =  0.7853981634F;  // π/4
static const float PI_OVER_3_f   =  1.0471975512F;  // π/3
static const float PI_OVER_2_f   =  1.57079632679F; // π/2
static const float TWO_PI_f      =  6.28318530718F; // 2π

static const char* SIMULATION_FILENAME = "simulation_data.json";
static const float SCALE = 1e7F; // scales in-sim meters to openGL coordinates -- this is an arbitrary number that can be adjusted... probably

// memory limits
#define MAX_PLANETS 32
#define PATH_CAPACITY 1024
#define MAX_SPACECRAFT 16
#define MAX_BURNS_PER_SPACECRAFT 16
#define MAX_NAME_LENGTH 32
#define JSON_BUFFER_SIZE 65536
#define SHADER_BUFFER_SIZE 16384
#define MAX_SPHERE_VERTICES 13500
#define MAX_LINE_BATCH 2048
#define MAX_FONT_CHARS 512
#define MAX_PAUSE_KEYFRAMES 1024

#define MAX_ERR_SIZE 128


#endif
