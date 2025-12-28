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

// sphere functions
sphere_mesh_t generateUnitSphere(unsigned int stacks, unsigned int sectors);
void freeSphere(sphere_mesh_t* sphere);

// text rendering functions
text_renderer_t initTextRenderer(const char* fontPath, unsigned int fontSize, int screenWidth, int screenHeight);
void renderText(text_renderer_t* renderer, const char* text, float x, float y, float scale, float color[3]);
void cleanupTextRenderer(text_renderer_t* renderer);

#endif //ORBITSIMULATION_GL_RENDERER_H