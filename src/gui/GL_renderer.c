//
// Created by java on 12/25/25.
//

#include "GL_renderer.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "../globals.h"
#include <ft2build.h>
#include FT_FREETYPE_H
#include "../math/matrix.h"

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

GLuint createShaderProgram(const char* vertexPath, const char* fragmentPath) {
    // load shader sources from files
    char* vertexShaderSource = loadShaderSource(vertexPath);
    char* fragmentShaderSource = loadShaderSource(fragmentPath);

    // compile the vertex shader
    GLuint vertexShader = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vertexShader, 1, (const char**)&vertexShaderSource, NULL);
    glCompileShader(vertexShader);

    // compile the fragment shader
    GLuint fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fragmentShader, 1, (const char**)&fragmentShaderSource, NULL);
    glCompileShader(fragmentShader);

    // link the shaders into a program
    GLuint shaderProgram = glCreateProgram();
    glAttachShader(shaderProgram, vertexShader);
    glAttachShader(shaderProgram, fragmentShader);
    glLinkProgram(shaderProgram);
    glDeleteShader(vertexShader);
    glDeleteShader(fragmentShader);

    // free the shader source strings (byebye)
    free(vertexShaderSource);
    free(fragmentShaderSource);

    return shaderProgram;
}

VBO_t createVBO(const float* vertices, const size_t vertexDataSize) {
    VBO_t vbo = {0};

    glGenVertexArrays(1, &vbo.VAO);
    glGenBuffers(1, &vbo.VBO);

    glBindVertexArray(vbo.VAO);
    glBindBuffer(GL_ARRAY_BUFFER, vbo.VBO);
    glBufferData(GL_ARRAY_BUFFER, (long)vertexDataSize, vertices, GL_STATIC_DRAW);

    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);

    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)(3 * sizeof(float)));
    glEnableVertexAttribArray(1);

    // unbind
    glBindVertexArray(0);

    return vbo;
}

void deleteVBO(const VBO_t vbo) {
    glDeleteVertexArrays(1, &vbo.VAO);
    glDeleteBuffers(1, &vbo.VBO);
}

// normalize a 3D vector
static void normalize_3d(float v[3]) {
    float length = sqrtf(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    if (length > 0.0f) {
        v[0] /= length;
        v[1] /= length;
        v[2] /= length;
    }
}

// cross product of two 3D vectors
static void cross_product_3d(const float a[3], const float b[3], float result[3]) {
    result[0] = a[1] * b[2] - a[2] * b[1];
    result[1] = a[2] * b[0] - a[0] * b[2];
    result[2] = a[0] * b[1] - a[1] * b[0];
}

// creates a view matrix that looks at the origin from cameraPos
mat4 createViewMatrix_originCentered(const float cameraPos[3]) {
    mat4 viewMatrix;
    float forward[3] = {-cameraPos[0], -cameraPos[1], -cameraPos[2]};
    normalize_3d(forward);
    float up[3] = {0.0f, 1.0f, 0.0f};
    float right[3];
    cross_product_3d(forward, up, right);
    normalize_3d(right);
    cross_product_3d(right, forward, up);

    // build view matrix
    viewMatrix.m[0] = right[0];
    viewMatrix.m[1] = up[0];
    viewMatrix.m[2] = -forward[0];
    viewMatrix.m[3] = 0.0f;
    viewMatrix.m[4] = right[1];
    viewMatrix.m[5] = up[1];
    viewMatrix.m[6] = -forward[1];
    viewMatrix.m[7] = 0.0f;
    viewMatrix.m[8] = right[2];
    viewMatrix.m[9] = up[2];
    viewMatrix.m[10] = -forward[2];
    viewMatrix.m[11] = 0.0f;
    viewMatrix.m[12] = -(right[0] * cameraPos[0] + right[1] * cameraPos[1] + right[2] * cameraPos[2]);
    viewMatrix.m[13] = -(up[0] * cameraPos[0] + up[1] * cameraPos[1] + up[2] * cameraPos[2]);
    viewMatrix.m[14] = (forward[0] * cameraPos[0] + forward[1] * cameraPos[1] + forward[2] * cameraPos[2]);
    viewMatrix.m[15] = 1.0f;
    return viewMatrix;
}

// creates a perspective projection matrix
mat4 createProjectionMatrix(const float fov, const float aspect, const float near, const float far) {
    mat4 projMatrix;
    float f = 1.0f / tanf(fov * 0.5f);

    projMatrix.m[0] = f / aspect;
    projMatrix.m[1] = 0.0f;
    projMatrix.m[2] = 0.0f;
    projMatrix.m[3] = 0.0f;

    projMatrix.m[4] = 0.0f;
    projMatrix.m[5] = f;
    projMatrix.m[6] = 0.0f;
    projMatrix.m[7] = 0.0f;

    projMatrix.m[8] = 0.0f;
    projMatrix.m[9] = 0.0f;
    projMatrix.m[10] = (far + near) / (near - far);
    projMatrix.m[11] = -1.0f;

    projMatrix.m[12] = 0.0f;
    projMatrix.m[13] = 0.0f;
    projMatrix.m[14] = (2.0f * far * near) / (near - far);
    projMatrix.m[15] = 0.0f;
    return projMatrix;
}

// sets a matrix uniform in the shader
void setMatrixUniform(const GLuint shaderProgram, const char* name, const mat4* matrix) {
    GLint location = glGetUniformLocation(shaderProgram, name);
    glUniformMatrix4fv(location, 1, GL_FALSE, matrix->m);
}

// main function for handling the camera, should be called in main().
void castCamera(const sim_properties_t sim, const GLuint shaderProgram) {
    // apply zoom to camera position unit vector
    float zoomedCameraPos[3] = {
        sim.wp.camera_pos.x * sim.wp.zoom,
        sim.wp.camera_pos.y * sim.wp.zoom,
        sim.wp.camera_pos.z * sim.wp.zoom
    };

    // create view matrix
    mat4 viewMatrix = createViewMatrix_originCentered(zoomedCameraPos);

    // create projection matrix
    float aspect = sim.wp.window_size_x / sim.wp.window_size_y;
    mat4 projMatrix = createProjectionMatrix(3.14159f / 4.0f, aspect, 0.1f, 100000.0f);

    // set matrices in shader
    setMatrixUniform(shaderProgram, "view", &viewMatrix);
    setMatrixUniform(shaderProgram, "projection", &projMatrix);
}

// sphere mesh generation
sphere_mesh_t generateUnitSphere(unsigned int stacks, unsigned int sectors) {
    sphere_mesh_t sphere = {0};

    // each quad (stack-sector intersection) creates 2 triangles = 6 vertices
    // total quads = stacks * sectors
    sphere.vertex_count = stacks * sectors * 6;
    sphere.data_size = sphere.vertex_count * 6 * sizeof(float); // 6 floats per vertex (pos + normal)

    sphere.vertices = (float*)malloc(sphere.data_size);
    if (!sphere.vertices) {
        fprintf(stderr, "Failed to allocate memory for sphere vertices\n");
        return sphere;
    }

    float* data = sphere.vertices;
    const float PI = 3.14159265359f;

    // generate vertices for each stack and sector
    for (unsigned int i = 0; i < stacks; ++i) {
        float theta1 = (float)i * PI / (float)stacks;
        float theta2 = (float)(i + 1) * PI / (float)stacks;

        for (unsigned int j = 0; j < sectors; ++j) {
            float phi1 = (float)j * 2.0f * PI / (float)sectors;
            float phi2 = (float)(j + 1) * 2.0f * PI / (float)sectors;

            // calculate 4 vertices of the quad
            // v1 (top-left)
            float v1_x = cosf(phi1) * sinf(theta1);
            float v1_y = cosf(theta1);
            float v1_z = sinf(phi1) * sinf(theta1);

            // v2 (bottom-left)
            float v2_x = cosf(phi1) * sinf(theta2);
            float v2_y = cosf(theta2);
            float v2_z = sinf(phi1) * sinf(theta2);

            // v3 (bottom-right)
            float v3_x = cosf(phi2) * sinf(theta2);
            float v3_y = cosf(theta2);
            float v3_z = sinf(phi2) * sinf(theta2);

            // v4 (top-right)
            float v4_x = cosf(phi2) * sinf(theta1);
            float v4_y = cosf(theta1);
            float v4_z = sinf(phi2) * sinf(theta1);

            // first triangle (v1, v2, v3)
            // v1
            *data++ = v1_x; *data++ = v1_y; *data++ = v1_z; // position
            *data++ = v1_x; *data++ = v1_y; *data++ = v1_z; // normal (same as position for unit sphere)
            // v2
            *data++ = v2_x; *data++ = v2_y; *data++ = v2_z;
            *data++ = v2_x; *data++ = v2_y; *data++ = v2_z;
            // v3
            *data++ = v3_x; *data++ = v3_y; *data++ = v3_z;
            *data++ = v3_x; *data++ = v3_y; *data++ = v3_z;

            // second triangle (v1, v3, v4)
            // v1
            *data++ = v1_x; *data++ = v1_y; *data++ = v1_z;
            *data++ = v1_x; *data++ = v1_y; *data++ = v1_z;
            // v3
            *data++ = v3_x; *data++ = v3_y; *data++ = v3_z;
            *data++ = v3_x; *data++ = v3_y; *data++ = v3_z;
            // v4
            *data++ = v4_x; *data++ = v4_y; *data++ = v4_z;
            *data++ = v4_x; *data++ = v4_y; *data++ = v4_z;
        }
    }

    return sphere;
}

void freeSphere(sphere_mesh_t* sphere) {
    if (sphere && sphere->vertices) {
        free(sphere->vertices);
        sphere->vertices = NULL;
        sphere->vertex_count = 0;
        sphere->data_size = 0;
    }
}

// line rendering stuff
line_batch_t createLineBatch(size_t max_lines) {
    line_batch_t batch = {0};

    batch.capacity = max_lines;
    batch.count = 0;

    // allocate vertex buffer
    size_t vertex_array_size = max_lines * 2 * 6 * sizeof(float);
    batch.vertices = (float*)malloc(vertex_array_size);

    if (!batch.vertices) {
        fprintf(stderr, "Failed to allocate memory for line batch\n");
        return batch;
    }

    // create VBO
    glGenVertexArrays(1, &batch.vbo.VAO);
    glGenBuffers(1, &batch.vbo.VBO);

    glBindVertexArray(batch.vbo.VAO);
    glBindBuffer(GL_ARRAY_BUFFER, batch.vbo.VBO);
    glBufferData(GL_ARRAY_BUFFER, (long)vertex_array_size, NULL, GL_DYNAMIC_DRAW);

    // position
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);

    // color
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)(3 * sizeof(float)));
    glEnableVertexAttribArray(1);

    glBindVertexArray(0);

    return batch;
}

void addLine(line_batch_t* batch, float x1, float y1, float z1, float x2, float y2, float z2, float r, float g, float b) {
    if (!batch || !batch->vertices || batch->count >= batch->capacity) {
        return;
    }

    // calculate index for this line's data
    size_t index = batch->count * 12; // 12 floats per line (2 vertices * 6 floats)

    // first vertex
    batch->vertices[index + 0] = x1;
    batch->vertices[index + 1] = y1;
    batch->vertices[index + 2] = z1;
    batch->vertices[index + 3] = r;
    batch->vertices[index + 4] = g;
    batch->vertices[index + 5] = b;

    // second vertex
    batch->vertices[index + 6] = x2;
    batch->vertices[index + 7] = y2;
    batch->vertices[index + 8] = z2;
    batch->vertices[index + 9] = r;
    batch->vertices[index + 10] = g;
    batch->vertices[index + 11] = b;

    batch->count++;
}

void renderLines(line_batch_t* batch, GLuint shader_program) {
    if (!batch || !batch->vertices || batch->count == 0) {
        return;
    }

    // activate the shader program
    glUseProgram(shader_program);

    // update VBO with all line data
    glBindBuffer(GL_ARRAY_BUFFER, batch->vbo.VBO);
    glBufferSubData(GL_ARRAY_BUFFER, 0, (long)(batch->count * 12 * sizeof(float)), batch->vertices);

    glBindVertexArray(batch->vbo.VAO);

    mat4 identity_mat = mat4_identity();
    setMatrixUniform(shader_program, "model", &identity_mat);

    // draw all lines in one call
    glDrawArrays(GL_LINES, 0, (GLsizei)(batch->count * 2));

    // reset for next frame
    batch->count = 0;
}

void freeLines(line_batch_t* batch) {
    if (!batch) {
        return;
    }

    if (batch->vertices) {
        free(batch->vertices);
        batch->vertices = NULL;
    }

    deleteVBO(batch->vbo);

    batch->capacity = 0;
    batch->count = 0;
}

// text rendering implementation
text_renderer_t initTextRenderer(const char* fontPath, unsigned int fontSize, int screenWidth, int screenHeight) {
    text_renderer_t renderer = {0};

    // initialize FreeType library
    FT_Library ft;
    if (FT_Init_FreeType(&ft)) {
        fprintf(stderr, "ERROR::FREETYPE: Could not init FreeType Library\n");
        return renderer;
    }

    // load font
    FT_Face face;
    if (FT_New_Face(ft, fontPath, 0, &face)) {
        fprintf(stderr, "ERROR::FREETYPE: Failed to load font from %s\n", fontPath);
        FT_Done_FreeType(ft);
        return renderer;
    }

    // set font size
    FT_Set_Pixel_Sizes(face, 0, fontSize);

    // disable byte-alignment restriction in OpenGL
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);

    // load first 128 ASCII characters
    for (unsigned char c = 0; c < 128; c++) {
        // load character glyph
        if (FT_Load_Char(face, c, FT_LOAD_RENDER)) {
            fprintf(stderr, "ERROR::FREETYPE: Failed to load Glyph for character %c\n", c);
            continue;
        }

        // generate texture
        GLuint texture;
        glGenTextures(1, &texture);
        glBindTexture(GL_TEXTURE_2D, texture);
        glTexImage2D(
            GL_TEXTURE_2D,
            0,
            GL_RED,
            (GLint)face->glyph->bitmap.width,
            (GLint)face->glyph->bitmap.rows,
            0,
            GL_RED,
            GL_UNSIGNED_BYTE,
            face->glyph->bitmap.buffer
        );

        // set texture options
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

        // store character
        character_t character = {
            texture,
            (GLint)face->glyph->bitmap.width,
            (GLint)face->glyph->bitmap.rows,
            face->glyph->bitmap_left,
            face->glyph->bitmap_top,
            (unsigned int)face->glyph->advance.x
        };
        renderer.characters[c] = character;
    }

    glBindTexture(GL_TEXTURE_2D, 0);

    // cleanup FreeType resources
    FT_Done_Face(face);
    FT_Done_FreeType(ft);

    // configure VAO/VBO for texture quads
    glGenVertexArrays(1, &renderer.VAO);
    glGenBuffers(1, &renderer.VBO);
    glBindVertexArray(renderer.VAO);
    glBindBuffer(GL_ARRAY_BUFFER, renderer.VBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(float) * 6 * 4, NULL, GL_DYNAMIC_DRAW);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 4, GL_FLOAT, GL_FALSE, 4 * sizeof(float), 0);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindVertexArray(0);

    // create shader program for text rendering
    const char* vertexShaderSource =
        "#version 330 core\n"
        "layout (location = 0) in vec4 vertex; // <vec2 pos, vec2 tex>\n"
        "out vec2 TexCoords;\n"
        "uniform mat4 projection;\n"
        "void main() {\n"
        "    gl_Position = projection * vec4(vertex.xy, 0.0, 1.0);\n"
        "    TexCoords = vertex.zw;\n"
        "}\0";

    const char* fragmentShaderSource =
        "#version 330 core\n"
        "in vec2 TexCoords;\n"
        "out vec4 color;\n"
        "uniform sampler2D text;\n"
        "uniform vec3 textColor;\n"
        "void main() {\n"
        "    vec4 sampled = vec4(1.0, 1.0, 1.0, texture(text, TexCoords).r);\n"
        "    color = vec4(textColor, 1.0) * sampled;\n"
        "}\0";

    // compile vertex shader
    GLuint vertexShader = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vertexShader, 1, &vertexShaderSource, NULL);
    glCompileShader(vertexShader);

    // compile fragment shader
    GLuint fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fragmentShader, 1, &fragmentShaderSource, NULL);
    glCompileShader(fragmentShader);

    // link shaders
    renderer.shader_program = glCreateProgram();
    glAttachShader(renderer.shader_program, vertexShader);
    glAttachShader(renderer.shader_program, fragmentShader);
    glLinkProgram(renderer.shader_program);

    glDeleteShader(vertexShader);
    glDeleteShader(fragmentShader);

    // set up orthographic projection matrix
    glUseProgram(renderer.shader_program);
    float projection[16] = {
        2.0f / (float)screenWidth, 0.0f, 0.0f, 0.0f,
        0.0f, 2.0f / (float)screenHeight, 0.0f, 0.0f,
        0.0f, 0.0f, -1.0f, 0.0f,
        -1.0f, -1.0f, 0.0f, 1.0f
    };
    GLint projLoc = glGetUniformLocation(renderer.shader_program, "projection");
    glUniformMatrix4fv(projLoc, 1, GL_FALSE, projection);

    return renderer;
}

void renderText(text_renderer_t* renderer, const char* text, float x, float y, float scale, float color[3]) {
    if (!renderer || !text) {
        return;
    }

    // enable blending for text transparency
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    // activate shader and set color
    glUseProgram(renderer->shader_program);
    glUniform3f(glGetUniformLocation(renderer->shader_program, "textColor"), color[0], color[1], color[2]);
    glActiveTexture(GL_TEXTURE0);
    glBindVertexArray(renderer->VAO);

    // iterate through all characters
    const char* c_ptr = text;
    while (*c_ptr) {
        character_t ch = renderer->characters[(unsigned char)*c_ptr];

        float xpos = x + (float)ch.bearing_x * scale;
        float ypos = y - (float)(ch.height - ch.bearing_y) * scale;

        float w = (float)ch.width * scale;
        float h = (float)ch.height * scale;

        // update VBO for each character
        float vertices[6][4] = {
            { xpos,     ypos + h,   0.0f, 0.0f },
            { xpos,     ypos,       0.0f, 1.0f },
            { xpos + w, ypos,       1.0f, 1.0f },

            { xpos,     ypos + h,   0.0f, 0.0f },
            { xpos + w, ypos,       1.0f, 1.0f },
            { xpos + w, ypos + h,   1.0f, 0.0f }
        };

        // render glyph texture over quad
        glBindTexture(GL_TEXTURE_2D, ch.texture_id);

        // update content of VBO memory
        glBindBuffer(GL_ARRAY_BUFFER, renderer->VBO);
        glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(vertices), vertices);
        glBindBuffer(GL_ARRAY_BUFFER, 0);

        // render quad
        glDrawArrays(GL_TRIANGLES, 0, 6);

        // advance cursor for next glyph (advance is in 1/64th pixels)
        x += (float)(ch.advance >> 6) * scale;
        c_ptr++;
    }

    glBindVertexArray(0);
    glBindTexture(GL_TEXTURE_2D, 0);
    glDisable(GL_BLEND);
}

void cleanupTextRenderer(text_renderer_t* renderer) {
    if (!renderer) {
        return;
    }

    // delete all character textures
    for (int i = 0; i < 128; i++) {
        if (renderer->characters[i].texture_id != 0) {
            glDeleteTextures(1, &renderer->characters[i].texture_id);
        }
    }

    // delete VAO and VBO
    if (renderer->VAO != 0) {
        glDeleteVertexArrays(1, &renderer->VAO);
    }
    if (renderer->VBO != 0) {
        glDeleteBuffers(1, &renderer->VBO);
    }

    // delete shader program
    if (renderer->shader_program != 0) {
        glDeleteProgram(renderer->shader_program);
    }
}

// render the coordinate plane to the screen
void renderCoordinatePlane(sim_properties_t sim, line_batch_t* line_batch) {
    float scale = sim.wp.zoom;

    // X axis (red)
    addLine(line_batch, -10.0f * scale, 0.0f, 0.0f, 10.0f * scale, 0.0f, 0.0f, 0.3f, 0.0f, 0.0f);

    // Y axis (green)
    addLine(line_batch, 0.0f, -10.0f * scale, 0.0f, 0.0f, 10.0f * scale, 0.0f, 0.0f, 0.3f, 0.0f);

    // Z axis (blue)
    addLine(line_batch, 0.0f, 0.0f, -10.0f * scale, 0.0f, 0.0f, 10.0f * scale, 0.0f, 0.0f, 0.3f);

    // perspective lines (gray)
    addLine(line_batch, 10.0f * scale, 0.0f, -10.0f * scale, -10.0f * scale, 0.0f, 10.0f * scale, 0.3f, 0.3f, 0.3f);
    addLine(line_batch, 10.0f * scale, 0.0f, 10.0f * scale, -10.0f * scale, 0.0f, -10.0f * scale, 0.3f, 0.3f, 0.3f);
}

// render the sim planets to the screen
void renderPlanets(sim_properties_t sim, GLuint shader_program, VBO_t planet_shape_buffer) {
    glBindVertexArray(planet_shape_buffer.VAO);
    for (int i = 0; i < sim.gb.count; i++) {
        // create a scale matrix based on the radius of the planet (pulled from JSON data)
        float size_scale_factor = (float)sim.gb.radius[i] / SCALE;
        mat4 scale_mat = mat4_scale(size_scale_factor, size_scale_factor, size_scale_factor);
        // create a translation matrix based on the current in-sim-world position of the planet
        mat4 translate_mat = mat4_translation((float)sim.gb.pos_x[i] / SCALE, (float)sim.gb.pos_y[i] / SCALE, (float)sim.gb.pos_z[i] / SCALE);
        // multiply the two matrices together to get a final scale/position matrix for the planet model on the screen
        mat4 planet_model = mat4_mul(translate_mat, scale_mat);
        // apply matrix and render to screen
        setMatrixUniform(shader_program, "model", &planet_model);
        glDrawArrays(GL_TRIANGLES, 0, sim.wp.planet_model_vertex_count);
    }
}

// render the sim crafts to the renderer
void renderCrafts(sim_properties_t sim, GLuint shader_program, VBO_t craft_shape_buffer) {
    glBindVertexArray(craft_shape_buffer.VAO);
    for (int i = 0; i < sim.gs.count; i++) {
        // create a scale matrix
        float size_scale_factor = 0.1f; // arbitrary scale based on what looks nice on screen
        mat4 scale_mat = mat4_scale(size_scale_factor, size_scale_factor * 2, size_scale_factor);
        // create a translation matrix based on the current in-sim-world position of the spacecraft
        mat4 translate_mat = mat4_translation((float)sim.gs.pos_x[i] / SCALE, (float)sim.gs.pos_y[i] / SCALE, (float)sim.gs.pos_z[i] / SCALE);
        // multiply the two matrices together to get a final scale/position matrix for the planet model on the screen
        mat4 spacecraft_model = mat4_mul(translate_mat, scale_mat);
        // apply matrix and render to screen
        setMatrixUniform(shader_program, "model", &spacecraft_model);
        glDrawArrays(GL_TRIANGLES, 0, 48);
    }
}

// render the stats on the screen
void renderStats(sim_properties_t sim, text_renderer_t text_renderer) {

    float white_color[] = { 1.0f, 1.0f, 1.0f};

    // calculate proper line height
    float line_height = (float)text_renderer.characters['H'].height * 1.5f;
    float cursor_starting_pos[2] = { 10.0f, (float)sim.wp.window_size_y - line_height - 10.0f}; // starting pos of writing
    float cursor_pos[2] = { cursor_starting_pos[0], cursor_starting_pos[1] };

    // text buffer used to hold text written to the stats window
    char text_buffer[32];

    // write planet velocities
    for (int i = 0; i < sim.gb.count; i++) {
        snprintf(text_buffer, sizeof(text_buffer), "%s's velocity: %.3f", sim.gb.names[i], sim.gb.vel[i]);
        renderText(&text_renderer, text_buffer, cursor_pos[0], cursor_pos[1], 1.0f, white_color);
        cursor_pos[1] -= line_height;
    }

}