//
// Created by java on 12/25/25.
//

#ifndef ORBITSIMULATION_GL_RENDERER_H
#define ORBITSIMULATION_GL_RENDERER_H

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
mat4 createViewMatrix(vec3_f cameraPos, vec3_f target);
mat4 createProjectionMatrix(float fov, float aspect, float near_plane, float far_plane);
void GL_setMatrixUniform(GLuint shaderProgram, const char* name, const mat4* matrix);
void castCamera(sim_properties_t sim, GLuint shaderProgram);

// drawing functions
sphere_mesh_t generateUnitSphere(unsigned int stacks, unsigned int sectors);
void freeSphere(sphere_mesh_t* sphere);

// line rendering
line_batch_t createLineBatch(size_t max_lines);
void addLine(line_batch_t* batch, float start_x, float start_y, float start_z, float end_x, float end_y, float end_z, float red, float green, float blue);
void renderLines(line_batch_t* batch, GLuint shader_program);
void freeLines(line_batch_t* batch);

font_t initFont(const char* path, float size);
void addText(font_t* font, float pos_x, float pos_y, const char* text, float scale);
void renderText(font_t* font, float window_width, float window_height, float red, float green, float blue);
void freeFont(const font_t* font);


// features (should come last)
void renderCoordinatePlane(sim_properties_t sim, line_batch_t* line_batch);
void renderPlanets(sim_properties_t sim, GLuint shader_program, VBO_t planet_shape_buffer);
void renderCrafts(sim_properties_t sim, GLuint shader_program, VBO_t craft_shape_buffer);
void renderStats(sim_properties_t sim, font_t* font);
void renderVisuals(sim_properties_t sim, line_batch_t* line_batch, craft_path_storage_t* craft_paths);

#endif //ORBITSIMULATION_GL_RENDERER_H