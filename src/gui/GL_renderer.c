//
// Created by java on 12/25/25.
//

#include "GL_renderer.h"
#include <stdio.h>
#include <stdlib.h>

char* loadShaderSource(const char* filepath) {
    FILE* file = fopen(filepath, "r");
    if (!file) {
        fprintf(stderr, "Failed to open shader file: %s\n", filepath);
        return NULL;
    }

    fseek(file, 0, SEEK_END);
    long size = ftell(file);
    fseek(file, 0, SEEK_SET);

    char* source = (char*)malloc(size + 1);

    fread(source, 1, size, file);
    source[size] = '\0';

    fclose(file);
    return source;
}

void worldToNDC(double world_x, double world_y,
                double camera_x, double camera_y,
                double meters_per_pixel,
                float screen_width, float screen_height,
                float* out_x, float* out_y) {
    // 1: transform world coordinates (meters) to camera-relative coordinates
    double rel_x = world_x - camera_x;
    double rel_y = world_y - camera_y;

    // 2: convert to pixel coordinates
    double pixel_x = rel_x / meters_per_pixel;
    double pixel_y = rel_y / meters_per_pixel;

    // 3: convert pixel coordinates to the magic openGL world
    // In opengl, (0,0) is center, (-1,-1) is bottom-left, (1,1) is top-right
    *out_x = (float)(2.0 * pixel_x / screen_width);
    *out_y = (float)(2.0 * pixel_y / screen_height);
}