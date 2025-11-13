#ifndef CALCULATION_FUNCTIONS_H
#define CALCULATION_FUNCTIONS_H

#include <stdio.h>
#include <math.h>
#include <SDL3/SDL.h>
#include <SDL3/SDL_main.h>
#include <SDL3_ttf/SDL_ttf.h>
#include <string.h>
#include <stdbool.h>
#include <stdlib.h>


extern const double G;
extern const int WINDOW_SIZE_X;
extern const int WINDOW_SIZE_Y;
extern const int ORIGIN_X;
extern const int ORIGIN_Y;
extern double meters_per_pixel;
extern const int FONT_SIZE;
extern TTF_Font* g_font;

typedef struct {
    double mass;
    int radius; // this is an approximate value just for visual purposes
    int visual_radius;
    double r_from_body;
    double pos_x;
    double pos_y;
    double vel_x;
    double vel_y;
    double vel;
    double acc_x;
    double acc_y;
    double force_x;
    double force_y;
    int pixel_coordinates_x;
    int pixel_coordinates_y;
} body_properties_t;

extern body_properties_t *global_bodies;
extern int num_bodies;

typedef struct {
    int x, y, width, height;
    bool is_hovered;
} speed_control_t;

void calculateForce(body_properties_t *b, body_properties_t b2);
void updateMotion(body_properties_t *b, double dt);
void transformCoordinates(body_properties_t *b);
void SDL_RenderFillCircle(SDL_Renderer* renderer, int centerX, int centerY, int radius);
void drawScaleBar(SDL_Renderer* renderer, double meters_per_pixel, int window_width, int window_height);
int calculateVisualRadius(body_properties_t body);
void SDL_WriteText(SDL_Renderer* renderer, TTF_Font* font, const char* text, float x, float y, SDL_Color color);

bool isMouseInRect(int mouse_x, int mouse_y, int rect_x, int rect_y, int rect_w, int rect_h);
void drawSpeedControl(SDL_Renderer* renderer, speed_control_t* control, double multiplier);
void runEventCheck(SDL_Event* event, bool* loop_running_condition, speed_control_t* speed_control, double* TIME_STEP, double* meters_per_pixel);

void drawStatsBox(SDL_Renderer* renderer, body_properties_t* bodies, int num_bodies, double sim_time);

void addOrbitalBody(double mass, double x_pos, double y_pos, double x_vel, double y_vel);

#endif