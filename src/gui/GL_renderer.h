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
mat4 createProjectionMatrix(const float fov, const float aspect, const float near_plane, const float far_plane);
void GL_setMatrixUniform(GLuint shaderProgram, const char* name, const mat4* matrix);
void castCamera(sim_properties_t sim, GLuint shaderProgram);

// drawing functions
sphere_mesh_t generateUnitSphere(unsigned int stacks, unsigned int sectors);
void freeSphere(sphere_mesh_t* sphere);

// line rendering
line_batch_t createLineBatch(size_t max_lines);
void addLine(line_batch_t* batch, const float start_x, const float start_y, const float start_z, const float end_x, const float end_y, const float end_z, const float red, const float green, const float blue);
void renderLines(line_batch_t* batch, GLuint shader_program);
void freeLines(line_batch_t* batch);

font_t initFont(const char* path, float size);
void addText(font_t* font, float pos_x, const float pos_y, const char* text, const float scale);
void renderText(font_t* font, const float window_width, const float window_height, const float red, const float green, const float blue);
void freeFont(const font_t* font);


// features (should come last)
void renderCoordinatePlane(sim_properties_t sim, line_batch_t* line_batch);
void renderPlanets(sim_properties_t sim, GLuint shader_program, VBO_t planet_shape_buffer);
void renderCrafts(sim_properties_t sim, GLuint shader_program, VBO_t craft_shape_buffer);
void renderStats(sim_properties_t sim, font_t* font);
void renderVisuals(sim_properties_t sim, line_batch_t* line_batch, craft_path_storage_t* craft_paths);

#endif //ORBITSIMULATION_GL_RENDERER_H