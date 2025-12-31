//
// Created by java on 12/25/25.
//

#ifndef ORBITSIMULATION_GL_RENDERER_H
#define ORBITSIMULATION_GL_RENDERER_H

#include <SDL3/SDL.h>
#include <GL/glew.h>
#ifdef __APPLE__
#include <OpenGL/gl.h>
#else
#include <GL/gl.h>
#endif
#include "../types.h"

char* loadShaderSource(const char* filepath);
GLuint createShaderProgram(const char* vertexPath, const char* fragmentPath);
VBO_t createVBO(const float* vertices, size_t vertexDataSize);
void deleteVBO(VBO_t vbo);

// camera helper functions
mat4 createViewMatrix_originCentered(const float cameraPos[3]);
mat4 createProjectionMatrix(float fov, float aspect, float near, float far);
void setMatrixUniform(GLuint shaderProgram, const char* name, const mat4* matrix);
void castCamera(sim_properties_t sim, GLuint shaderProgram);

// drawing functions
sphere_mesh_t generateUnitSphere(unsigned int stacks, unsigned int sectors);
void freeSphere(sphere_mesh_t* sphere);

// line rendering
line_batch_t createLineBatch(size_t max_lines);
void addLine(line_batch_t* batch, float x1, float y1, float z1, float x2, float y2, float z2, float r, float g, float b);
void renderLines(line_batch_t* batch, GLuint shader_program);
void freeLines(line_batch_t* batch);

// text rendering
typedef struct {
    GLuint tex, shader, vao, vbo;
    float* verts;
    int count;
} font_t;

font_t initFont(const char* path, float size);
void addText(font_t* font, float x, float y, const char* text, float scale);
void renderText(font_t* font, float window_w, float window_h, float r, float g, float b);
void freeFont(font_t* font);


// features (should come last)
void renderCoordinatePlane(sim_properties_t sim, line_batch_t* line_batch);
void renderPlanets(sim_properties_t sim, GLuint shader_program, VBO_t planet_shape_buffer);
void renderCrafts(sim_properties_t sim, GLuint shader_program, VBO_t craft_shape_buffer);
void renderStats(sim_properties_t sim, font_t* font);
void renderDebug(sim_properties_t sim, line_batch_t* line_batch);

#endif //ORBITSIMULATION_GL_RENDERER_H