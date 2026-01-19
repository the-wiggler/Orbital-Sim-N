//
// Created by java on 12/25/25.
//

#include "GL_renderer.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../globals.h"
#include "../math/matrix.h"

char* loadShaderSource(const char* filepath) {
    FILE* file = fopen(filepath, "rb");
    if (!file) {
        fprintf(stderr, "Failed to open shader file: %s\n", filepath);
        return NULL;
    }

    fseek(file, 0, SEEK_END);
    const long size = ftell(file);
    fseek(file, 0, SEEK_SET);

    char* source = (char*)malloc(size + 1);
    if (!source) {
        fprintf(stderr, "Failed to allocate memory for shader source\n");
        fclose(file);
        return NULL;
    }

    const size_t bytesRead = fread(source, 1, size, file);
    source[bytesRead] = '\0';  // Use actual bytes read, not file size

    fclose(file);
    return source;
}

GLuint createShaderProgram(const char* vertexPath, const char* fragmentPath) {
    // load shader sources from files
    char* vertexShaderSource = loadShaderSource(vertexPath);
    char* fragmentShaderSource = loadShaderSource(fragmentPath);

    if (!vertexShaderSource) {
        fprintf(stderr, "Failed to load vertex shader: %s\n", vertexPath);
        return 0;
    }
    if (!fragmentShaderSource) {
        fprintf(stderr, "Failed to load fragment shader: %s\n", fragmentPath);
        free(vertexShaderSource);
        return 0;
    }

    // compile the vertex shader
    const GLuint vertexShader = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vertexShader, 1, (const char**)&vertexShaderSource, NULL);
    glCompileShader(vertexShader);

    // check vertex shader compilation
    GLint success;
    GLchar infoLog[512];
    glGetShaderiv(vertexShader, GL_COMPILE_STATUS, &success);
    if (!success) {
        glGetShaderInfoLog(vertexShader, 512, NULL, infoLog);
        fprintf(stderr, "ERROR: Vertex shader compilation failed (%s):\n%s\n", vertexPath, infoLog);
        free(vertexShaderSource);
        free(fragmentShaderSource);
        return 0;
    }

    // compile the fragment shader
    const GLuint fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fragmentShader, 1, (const char**)&fragmentShaderSource, NULL);
    glCompileShader(fragmentShader);

    // check fragment shader compilation
    glGetShaderiv(fragmentShader, GL_COMPILE_STATUS, &success);
    if (!success) {
        glGetShaderInfoLog(fragmentShader, 512, NULL, infoLog);
        fprintf(stderr, "ERROR: Fragment shader compilation failed (%s):\n%s\n", fragmentPath, infoLog);
        glDeleteShader(vertexShader);
        free(vertexShaderSource);
        free(fragmentShaderSource);
        return 0;
    }

    // link the shaders into a program
    const GLuint shaderProgram = glCreateProgram();
    glAttachShader(shaderProgram, vertexShader);
    glAttachShader(shaderProgram, fragmentShader);
    glLinkProgram(shaderProgram);

    // check linking errors
    glGetProgramiv(shaderProgram, GL_LINK_STATUS, &success);
    if (!success) {
        glGetProgramInfoLog(shaderProgram, 512, NULL, infoLog);
        fprintf(stderr, "ERROR: Shader program linking failed:\n%s\n", infoLog);
        glDeleteShader(vertexShader);
        glDeleteShader(fragmentShader);
        free(vertexShaderSource);
        free(fragmentShaderSource);
        return 0;
    }

    glDeleteShader(vertexShader);
    glDeleteShader(fragmentShader);

    // free the shader source strings (byebye)
    free(vertexShaderSource);
    free(fragmentShaderSource);

    printf("Successfully compiled shader program: %s + %s\n", vertexPath, fragmentPath);
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

// creates a view matrix that looks at target from cameraPos
mat4 createViewMatrix(const vec3_f cameraPos, const vec3_f target) {
    mat4 viewMatrix;
    vec3_f forward = vec3_f_normalize(vec3_f_sub(target, cameraPos));
    vec3_f up = {0.0f, 0.0f, 1.0f}; // z is "UP"
    vec3_f right = vec3_f_normalize(vec3_f_cross(forward, up));
    up = vec3_f_cross(right, forward);

    // build view matrix
    viewMatrix.m[0] = right.x;
    viewMatrix.m[1] = up.x;
    viewMatrix.m[2] = -forward.x;
    viewMatrix.m[3] = 0.0f;
    viewMatrix.m[4] = right.y;
    viewMatrix.m[5] = up.y;
    viewMatrix.m[6] = -forward.y;
    viewMatrix.m[7] = 0.0f;
    viewMatrix.m[8] = right.z;
    viewMatrix.m[9] = up.z;
    viewMatrix.m[10] = -forward.z;
    viewMatrix.m[11] = 0.0f;
    viewMatrix.m[12] = -vec3_f_dot(right, cameraPos);
    viewMatrix.m[13] = -vec3_f_dot(up, cameraPos);
    viewMatrix.m[14] = vec3_f_dot(forward, cameraPos);
    viewMatrix.m[15] = 1.0f;
    return viewMatrix;
}

// creates a perspective projection matrix
mat4 createProjectionMatrix(const float fov, const float aspect, const float near, const float far) {
    mat4 projMatrix;
    const float f = 1.0f / tanf(fov * 0.5f);

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
void GL_setMatrixUniform(const GLuint shaderProgram, const char* name, const mat4* matrix) {
    const GLint location = glGetUniformLocation(shaderProgram, name);
    glUniformMatrix4fv(location, 1, GL_FALSE, matrix->m);
}

// main function for handling the camera, should be called in main().
void castCamera(const sim_properties_t sim, const GLuint shaderProgram) {
    // camera position = target + direction * zoom (orbits around target)
    const vec3_f cameraPos = {
        sim.wp.cam_target.x + sim.wp.camera_pos.x * sim.wp.zoom,
        sim.wp.cam_target.y + sim.wp.camera_pos.y * sim.wp.zoom,
        sim.wp.cam_target.z + sim.wp.camera_pos.z * sim.wp.zoom
    };

    // create view matrix
    const mat4 viewMatrix = createViewMatrix(cameraPos, sim.wp.cam_target);

    // create projection matrix
    const float aspect = sim.wp.window_size_x / sim.wp.window_size_y;
    const mat4 projMatrix = createProjectionMatrix(PI_OVER_4_f, aspect, 0.1f, 100000.0f);

    // set matrices in shader
    GL_setMatrixUniform(shaderProgram, "view", &viewMatrix);
    GL_setMatrixUniform(shaderProgram, "projection", &projMatrix);
}

// sphere mesh generation
sphere_mesh_t generateUnitSphere(const unsigned int stacks, const unsigned int sectors) {
    sphere_mesh_t sphere = {0};
    sphere.vertex_count = stacks * sectors * 6;
    sphere.data_size = sphere.vertex_count * 6 * sizeof(float);

    sphere.vertices = (float*)malloc(sphere.data_size);
    if (!sphere.vertices) {
        fprintf(stderr, "Failed to allocate memory for sphere vertices\n");
        return sphere;
    }

    float* data = sphere.vertices;

    // generate vertices for each stack and sector
    for (unsigned int i = 0; i < stacks; ++i) {
        const float theta1 = (float)i * PI_f / (float)stacks;
        const float theta2 = (float)(i + 1) * PI_f / (float)stacks;

        for (unsigned int j = 0; j < sectors; ++j) {
            const float phi1 = (float)j * 2.0f * PI_f / (float)sectors;
            const float phi2 = (float)(j + 1) * 2.0f * PI_f / (float)sectors;
            const float v1_x = cosf(phi1) * sinf(theta1);
            const float v1_y = sinf(phi1) * sinf(theta1);
            const float v1_z = cosf(theta1);
            const float v2_x = cosf(phi1) * sinf(theta2);
            const float v2_y = sinf(phi1) * sinf(theta2);
            const float v2_z = cosf(theta2);
            const float v3_x = cosf(phi2) * sinf(theta2);
            const float v3_y = sinf(phi2) * sinf(theta2);
            const float v3_z = cosf(theta2);
            const float v4_x = cosf(phi2) * sinf(theta1);
            const float v4_y = sinf(phi2) * sinf(theta1);
            const float v4_z = cosf(theta1);
            *data++ = v1_x; *data++ = v1_y; *data++ = v1_z;
            *data++ = v1_x; *data++ = v1_y; *data++ = v1_z;
            *data++ = v2_x; *data++ = v2_y; *data++ = v2_z;
            *data++ = v2_x; *data++ = v2_y; *data++ = v2_z;
            *data++ = v3_x; *data++ = v3_y; *data++ = v3_z;
            *data++ = v3_x; *data++ = v3_y; *data++ = v3_z;
            *data++ = v1_x; *data++ = v1_y; *data++ = v1_z;
            *data++ = v1_x; *data++ = v1_y; *data++ = v1_z;
            *data++ = v3_x; *data++ = v3_y; *data++ = v3_z;
            *data++ = v3_x; *data++ = v3_y; *data++ = v3_z;
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
line_batch_t createLineBatch(const size_t max_lines) {
    line_batch_t batch = {0};

    batch.capacity = max_lines;
    batch.count = 0;

    // allocate vertex buffer
    const size_t vertex_array_size = max_lines * 2 * 6 * sizeof(float);
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

void addLine(line_batch_t* batch, const float x1, const float y1, const float z1, const float x2, const float y2, const float z2, const float r, const float g, const float b) {
    if (!batch || !batch->vertices || batch->count >= batch->capacity) {
        return;
    }

    // calculate index for this line's data
    const size_t index = batch->count * 12; // 12 floats per line (2 vertices * 6 floats)

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

void renderLines(line_batch_t* batch, const GLuint shader_program) {
    if (!batch || !batch->vertices || batch->count == 0) {
        return;
    }

    // activate the shader program
    glUseProgram(shader_program);

    // update VBO with all line data
    glBindBuffer(GL_ARRAY_BUFFER, batch->vbo.VBO);
    glBufferSubData(GL_ARRAY_BUFFER, 0, (long)(batch->count * 12 * sizeof(float)), batch->vertices);

    glBindVertexArray(batch->vbo.VAO);

    const mat4 identity_mat = mat4_identity();
    GL_setMatrixUniform(shader_program, "model", &identity_mat);

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

#define STB_TRUETYPE_IMPLEMENTATION
#include "stb_truetype.h"

#define ATLAS 512
#define MAX_CHARS 4096

static stbtt_bakedchar cdata[96];

font_t initFont(const char* path, const float size) {
    font_t f = {0};

    FILE* fp = fopen(path, "rb");
    if (!fp) {
        fprintf(stderr, "Failed to open font file: %s\n", path);
        return f;
    }
    fseek(fp, 0, SEEK_END);
    const long len = ftell(fp);
    rewind(fp);
    unsigned char* buf = malloc(len);
    fread(buf, 1, len, fp);
    fclose(fp);

    unsigned char* bmp = malloc(ATLAS * ATLAS);
    stbtt_BakeFontBitmap(buf, 0, size, bmp, ATLAS, ATLAS, 32, 96, cdata);
    free(buf);

    glGenTextures(1, &f.tex);
    glBindTexture(GL_TEXTURE_2D, f.tex);

    GLint internalFormat = GL_RED;

#ifdef __EMSCRIPTEN__
    // WebGL needs a different internal format
    internalFormat = GL_R8;
#endif

    glTexImage2D(GL_TEXTURE_2D, 0, internalFormat, ATLAS, ATLAS, 0, GL_RED, GL_UNSIGNED_BYTE, bmp);

    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    free(bmp);

    f.shader = createShaderProgram("shaders/text.vert", "shaders/text.frag");
    if (f.shader == 0) {
        fprintf(stderr, "Failed to create text shader program\n");
        glDeleteTextures(1, &f.tex);
        return f;
    }

    f.verts = malloc(MAX_CHARS * 24 * sizeof(float));

    glGenVertexArrays(1, &f.vao);
    glGenBuffers(1, &f.vbo);
    glBindVertexArray(f.vao);
    glBindBuffer(GL_ARRAY_BUFFER, f.vbo);
    glBufferData(GL_ARRAY_BUFFER, MAX_CHARS * 24 * sizeof(float), NULL, GL_DYNAMIC_DRAW);
    glVertexAttribPointer(0, 4, GL_FLOAT, GL_FALSE, 4 * sizeof(float), 0);
    glEnableVertexAttribArray(0);

    printf("Successfully initialized font: %s\n", path);
    return f;
}

void addText(font_t* font, float x, const float y, const char* text, const float scale) {
    for (; *text; text++) {
        const int c = *text - 32;
        if (c < 0 || c >= 96) continue;
        const stbtt_bakedchar* b = &cdata[c];

        float x0 = x + b->xoff * scale, y0 = y + b->yoff * scale;
        float x1 = x0 + (float)(b->x1 - b->x0) * scale, y1 = y0 + (float)(b->y1 - b->y0) * scale;
        float s0 = (float)b->x0 / (float)ATLAS, t0 = (float)b->y0 / (float)ATLAS;
        float s1 = (float)b->x1 / (float)ATLAS, t1 = (float)b->y1 / (float)ATLAS;

        float* v = &font->verts[font->count++ * 24];
        v[0]=x0; v[1]=y0; v[2]=s0; v[3]=t0;
        v[4]=x1; v[5]=y0; v[6]=s1; v[7]=t0;
        v[8]=x1; v[9]=y1; v[10]=s1; v[11]=t1;
        v[12]=x0; v[13]=y0; v[14]=s0; v[15]=t0;
        v[16]=x1; v[17]=y1; v[18]=s1; v[19]=t1;
        v[20]=x0; v[21]=y1; v[22]=s0; v[23]=t1;

        x += b->xadvance * scale;
    }
}

void renderText(font_t* font, const float window_w, const float window_h, const float r, const float g, const float b) {
    if (!font->count) return;

    glDisable(GL_DEPTH_TEST);
    glUseProgram(font->shader);

    const float proj[16] =   {2/window_w,0,0,0,
                        0,-2/window_h,0,0,
                        0,0,-1,0,
                        -1,1,0,1};
    glUniformMatrix4fv(glGetUniformLocation(font->shader, "proj"), 1, GL_FALSE, proj);
    glUniform3f(glGetUniformLocation(font->shader, "color"), r, g, b);

    glBindTexture(GL_TEXTURE_2D, font->tex);
    glBindBuffer(GL_ARRAY_BUFFER, font->vbo);
    glBufferSubData(GL_ARRAY_BUFFER, 0, font->count * 24 * sizeof(float), font->verts);
    glBindVertexArray(font->vao);
    glDrawArrays(GL_TRIANGLES, 0, font->count * 6);

    glEnable(GL_DEPTH_TEST);
    font->count = 0;
}

void freeFont(const font_t* font) {
    free(font->verts);
    glDeleteTextures(1, &font->tex);
    glDeleteProgram(font->shader);
    glDeleteVertexArrays(1, &font->vao);
    glDeleteBuffers(1, &font->vbo);
}








// render the coordinate plane to the screen
void renderCoordinatePlane(const sim_properties_t sim, line_batch_t* line_batch) {
    const float scale = sim.wp.zoom;

    // X axis (red) - horizontal
    addLine(line_batch, -10.0f * scale, 0.0f, 0.0f, 10.0f * scale, 0.0f, 0.0f, 0.3f, 0.0f, 0.0f);

    // Y axis (green) - horizontal
    addLine(line_batch, 0.0f, -10.0f * scale, 0.0f, 0.0f, 10.0f * scale, 0.0f, 0.0f, 0.3f, 0.0f);

    // Z axis (blue) - vertical "up"
    addLine(line_batch, 0.0f, 0.0f, -10.0f * scale, 0.0f, 0.0f, 10.0f * scale, 0.0f, 0.0f, 0.3f);

    // perspective lines in X-Y plane (gray)
    addLine(line_batch, 10.0f * scale, 10.0f * scale, 0.0f, -10.0f * scale, -10.0f * scale, 0.0f, 0.3f, 0.3f, 0.3f);
    addLine(line_batch, 10.0f * scale, -10.0f * scale, 0.0f, -10.0f * scale, 10.0f * scale, 0.0f, 0.3f, 0.3f, 0.3f);
}

// render the sim planets to the screen
void renderPlanets(const sim_properties_t sim, const GLuint shader_program, const VBO_t planet_shape_buffer) {
    glBindVertexArray(planet_shape_buffer.VAO);
    for (int i = 0; i < sim.gb.count; i++) {
        const body_t* body = &sim.gb.bodies[i];

        // create a scale matrix based on the radius of the planet
        const float size_scale_factor = (float)body->radius / SCALE;
        const mat4 scale_mat = mat4_scale(size_scale_factor, size_scale_factor, size_scale_factor);

        // create a rotation matrix based on the planet's attitude quaternion
        const mat4 rotation_mat = quaternionToMatrix(body->attitude);

        // create a translation matrix based on the current position
        const mat4 translate_mat = mat4_translation(
            (float)body->pos.x / SCALE,
            (float)body->pos.y / SCALE,
            (float)body->pos.z / SCALE);

        // multiply matrices together
        const mat4 temp = mat4_mul(rotation_mat, scale_mat);
        mat4 planet_model = mat4_mul(translate_mat, temp);

        // apply matrix and render to screen
        GL_setMatrixUniform(shader_program, "model", &planet_model);
        glDrawArrays(GL_TRIANGLES, 0, sim.wp.planet_model_vertex_count);
    }

    // draw sphere of influence
    if (sim.wp.draw_planet_SOI) {
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        glDepthMask(GL_FALSE);
        glUniform1i(glGetUniformLocation(shader_program, "useOverride"), 1);
        glUniform4f(glGetUniformLocation(shader_program, "colorOverride"), 0.0f, 0.5f, 1.0f, 0.1f);
        glBindVertexArray(planet_shape_buffer.VAO);
        for (int i = 0; i < sim.gb.count; i++) {
            const body_t* body = &sim.gb.bodies[i];
            const float size_scale_factor = (float)body->SOI_radius / SCALE;
            const mat4 scale_mat = mat4_scale(size_scale_factor, size_scale_factor, size_scale_factor);
            const mat4 rotation_mat = quaternionToMatrix(body->attitude);
            const mat4 translate_mat = mat4_translation((float)body->pos.x / SCALE,(float)body->pos.y / SCALE,(float)body->pos.z / SCALE);
            const mat4 temp = mat4_mul(rotation_mat, scale_mat);
            mat4 planet_model = mat4_mul(translate_mat, temp);
            GL_setMatrixUniform(shader_program, "model", &planet_model);
            glDrawArrays(GL_TRIANGLES, 0, sim.wp.planet_model_vertex_count);
        }
        glUniform1i(glGetUniformLocation(shader_program, "useOverride"), 0);
        glDepthMask(GL_TRUE);
        glDisable(GL_BLEND);
    }
}

// render the sim crafts to the renderer
void renderCrafts(const sim_properties_t sim, const GLuint shader_program, const VBO_t craft_shape_buffer) {
    glBindVertexArray(craft_shape_buffer.VAO);
    for (int i = 0; i < sim.gs.count; i++) {
        const spacecraft_t* craft = &sim.gs.spacecraft[i];

        // create a scale matrix
        const float size_scale_factor = 0.05f;
        const mat4 S = mat4_scale(size_scale_factor, size_scale_factor * 2, size_scale_factor);
        const mat4 R = quaternionToMatrix(craft->attitude);
        const mat4 T = mat4_translation(
            (float)craft->pos.x / SCALE,
            (float)craft->pos.y / SCALE,
            (float)craft->pos.z / SCALE);

        const mat4 temp = mat4_mul(R, S);
        mat4 model = mat4_mul(T, temp);

        // apply matrix and render to screen
        GL_setMatrixUniform(shader_program, "model", &model);
        glDrawArrays(GL_TRIANGLES, 0, 48);
    }
}

// render the stats on the screen
void renderStats(const sim_properties_t sim, font_t* font) {

    // calculate proper line height
    const float line_height = 20.0f;
    const float cursor_starting_pos[2] = { 10.0f, line_height + 10.0f};
    float cursor_pos[2] = { cursor_starting_pos[0], cursor_starting_pos[1] };

    // text buffer used to hold text written to the stats window
    char text_buffer[64];

    // paused indication
    if (sim.wp.sim_running) snprintf(text_buffer, sizeof(text_buffer), "Sim running");
    else snprintf(text_buffer, sizeof(text_buffer), "Sim paused");
    addText(font, cursor_pos[0], cursor_pos[1], text_buffer, 0.8f);
    cursor_pos[1] += line_height;

    // write time step
    snprintf(text_buffer, sizeof(text_buffer), "Step: %.4g", sim.wp.time_step);
    addText(font, cursor_pos[0], cursor_pos[1], text_buffer, 0.8f);
    cursor_pos[1] += line_height;

    // time indication
    const double time = sim.wp.sim_time / 3600;
    if (time < 72.0) snprintf(text_buffer, sizeof(text_buffer), "Time: %.2f hrs", time);
    else if (time < 8766.0) snprintf(text_buffer, sizeof(text_buffer), "Time: %.2f days", time / 24);
    addText(font, cursor_pos[0], cursor_pos[1], text_buffer, 0.8f);
    cursor_pos[1] += line_height;

    // spacer
    cursor_pos[1] += line_height;

    for (int i = 0; i < sim.gs.count; i++) {
        const int closest_id = sim.gs.spacecraft[i].closest_planet_id;
        const int soi_id = sim.gs.spacecraft[i].SOI_planet_id;

        snprintf(text_buffer, sizeof(text_buffer), "%s", sim.gs.spacecraft[i].name);
        addText(font, cursor_pos[0], cursor_pos[1], text_buffer, 0.8f);
        cursor_pos[1] += line_height;

        snprintf(text_buffer, sizeof(text_buffer), "Closest Planet: %s", sim.gb.bodies[closest_id].name);
        addText(font, cursor_pos[0], cursor_pos[1], text_buffer, 0.7f);
        cursor_pos[1] += line_height;

        snprintf(text_buffer, sizeof(text_buffer), "Distance: %.2f km", sqrt(sim.gs.spacecraft[i].closest_r_squared) / 1000 - sim.gb.bodies[closest_id].radius / 1000);
        addText(font, cursor_pos[0], cursor_pos[1], text_buffer, 0.7f);
        cursor_pos[1] += line_height;

        snprintf(text_buffer, sizeof(text_buffer), "In SOI of: %s", sim.gb.bodies[soi_id].name);
        addText(font, cursor_pos[0], cursor_pos[1], text_buffer, 0.7f);
        cursor_pos[1] += line_height;

        snprintf(text_buffer, sizeof(text_buffer), "Semi Major Axis: %.4g", sim.gs.spacecraft[i].semi_major_axis);
        addText(font, cursor_pos[0], cursor_pos[1], text_buffer, 0.7f);
        cursor_pos[1] += line_height;
        snprintf(text_buffer, sizeof(text_buffer), "Eccentricity: %.6f", sim.gs.spacecraft[i].eccentricity);
        addText(font, cursor_pos[0], cursor_pos[1], text_buffer, 0.7f);
        cursor_pos[1] += line_height;
        snprintf(text_buffer, sizeof(text_buffer), "Inclination: %.4f", sim.gs.spacecraft[i].inclination);
        addText(font, cursor_pos[0], cursor_pos[1], text_buffer, 0.7f);
        cursor_pos[1] += line_height;
        snprintf(text_buffer, sizeof(text_buffer), "Ascending Node: %.4f", sim.gs.spacecraft[i].ascending_node);
        addText(font, cursor_pos[0], cursor_pos[1], text_buffer, 0.7f);
        cursor_pos[1] += line_height;
        snprintf(text_buffer, sizeof(text_buffer), "Arg of Periapsis: %.4f", sim.gs.spacecraft[i].arg_periapsis);
        addText(font, cursor_pos[0], cursor_pos[1], text_buffer, 0.7f);
        cursor_pos[1] += line_height;
        snprintf(text_buffer, sizeof(text_buffer), "True Anomaly: %.4f", sim.gs.spacecraft[i].true_anomaly);
        addText(font, cursor_pos[0], cursor_pos[1], text_buffer, 0.7f);
        cursor_pos[1] += line_height;

        cursor_pos[1] += line_height;
    }
}

void renderPlanetPaths(sim_properties_t* sim, line_batch_t* line_batch, object_path_storage_t* planet_paths) {
    // initialize or resize planet paths if needed
    if (sim->gb.count > 0 && planet_paths->num_objects != sim->gb.count) {
        free(planet_paths->positions);
        free(planet_paths->counts);
        planet_paths->num_objects = sim->gb.count;
        planet_paths->capacity = PATH_CAPACITY;
        planet_paths->positions = malloc(planet_paths->num_objects * planet_paths->capacity * sizeof(vec3));
        planet_paths->counts = calloc(planet_paths->num_objects, sizeof(int));
    }

    if (sim->wp.draw_planet_path && planet_paths->num_objects > 0 && planet_paths->positions != NULL && planet_paths->counts != NULL) {
        // draw orbital paths for all planets
        for (int p = 0; p < planet_paths->num_objects; p++) {
            const int base = p * planet_paths->capacity;
            for (int i = 1; i < planet_paths->counts[p]; i++) {
                const vec3 p1 = planet_paths->positions[base + i - 1];
                const vec3 p2 = planet_paths->positions[base + i];
                addLine(line_batch, (float)p1.x, (float)p1.y, (float)p1.z,
                                    (float)p2.x, (float)p2.y, (float)p2.z, 0.5f, 1.0f, 0.5f);
            }
        }
    }

    // record planet paths
    if (sim->wp.frame_counter % 5 == 0 && sim->gb.count > 0 && planet_paths->counts != NULL && planet_paths->positions != NULL) {
        for (int p = 0; p < sim->gb.count; p++) {
            const body_t* body = &sim->gb.bodies[p];
            const int idx = p * planet_paths->capacity + planet_paths->counts[p];
            if (planet_paths->counts[p] < planet_paths->capacity) {
                planet_paths->positions[idx].x = body->pos.x / SCALE;
                planet_paths->positions[idx].y = body->pos.y / SCALE;
                planet_paths->positions[idx].z = body->pos.z / SCALE;
                planet_paths->counts[p]++;
            } else {
                const int base = p * planet_paths->capacity;
                for (int i = 1; i < planet_paths->capacity; i++) {
                    planet_paths->positions[base + i - 1] = planet_paths->positions[base + i];
                }
                planet_paths->positions[base + planet_paths->capacity - 1].x = body->pos.x / SCALE;
                planet_paths->positions[base + planet_paths->capacity - 1].y = body->pos.y / SCALE;
                planet_paths->positions[base + planet_paths->capacity - 1].z = body->pos.z / SCALE;
            }
        }
    }
}

void renderCraftPaths(sim_properties_t* sim, line_batch_t* line_batch, object_path_storage_t* craft_paths) {
    window_params_t wp = sim->wp;
    spacecraft_properties_t gs = sim->gs;

    if (wp.draw_craft_path && craft_paths->num_objects > 0) {
        // draw orbital paths for all crafts
        for (int p = 0; p < craft_paths->num_objects; p++) {
            const int base = p * craft_paths->capacity;
            for (int i = 1; i < craft_paths->counts[p]; i++) {
                const vec3 p1 = craft_paths->positions[base + i - 1];
                const vec3 p2 = craft_paths->positions[base + i];
                addLine(line_batch, (float)p1.x, (float)p1.y, (float)p1.z,
                                    (float)p2.x, (float)p2.y, (float)p2.z, 1.0f, 1.0f, 0.5f);
            }
        }
    }

    // initialize or resize craft paths if needed
    if (gs.count > 0 && craft_paths->num_objects != gs.count) {
        free(craft_paths->positions);
        free(craft_paths->counts);
        craft_paths->num_objects = gs.count;
        craft_paths->capacity = PATH_CAPACITY;
        craft_paths->positions = malloc(craft_paths->num_objects * craft_paths->capacity * sizeof(vec3));
        craft_paths->counts = calloc(craft_paths->num_objects, sizeof(int));
    }

    // record craft paths
    if (gs.count > 0) {
        for (int p = 0; p < gs.count; p++) {
            const spacecraft_t* craft = &sim->gs.spacecraft[p];
            const int idx = p * craft_paths->capacity + craft_paths->counts[p];
            if (craft_paths->counts[p] < craft_paths->capacity) {
                craft_paths->positions[idx].x = craft->pos.x / SCALE;
                craft_paths->positions[idx].y = craft->pos.y / SCALE;
                craft_paths->positions[idx].z = craft->pos.z / SCALE;
                craft_paths->counts[p]++;
            } else {
                const int base = p * craft_paths->capacity;
                for (int i = 1; i < craft_paths->capacity; i++) {
                    craft_paths->positions[base + i - 1] = craft_paths->positions[base + i];
                }
                craft_paths->positions[base + craft_paths->capacity - 1].x = craft->pos.x / SCALE;
                craft_paths->positions[base + craft_paths->capacity - 1].y = craft->pos.y / SCALE;
                craft_paths->positions[base + craft_paths->capacity - 1].z = craft->pos.z / SCALE;
            }
        }
    }
}

// renders debug features when they are enabled
void renderVisuals(sim_properties_t sim, line_batch_t* line_batch, object_path_storage_t* planet_paths, object_path_storage_t* craft_paths) {
    window_params_t wp = sim.wp;
    spacecraft_properties_t gs = sim.gs;
    body_properties_t gb = sim.gb;

    // create temporary arrays for scaled positions (avoids modifying original data)
    vec3_f scaled_body_pos[gb.count];
    vec3_f scaled_craft_pos[gs.count];

    for (int i = 0; i < gb.count; i++) {
        scaled_body_pos[i].x = (float)(gb.bodies[i].pos.x / SCALE);
        scaled_body_pos[i].y = (float)(gb.bodies[i].pos.y / SCALE);
        scaled_body_pos[i].z = (float)(gb.bodies[i].pos.z / SCALE);
    }
    for (int i = 0; i < gs.count; i++) {
        scaled_craft_pos[i].x = (float)(gs.spacecraft[i].pos.x / SCALE);
        scaled_craft_pos[i].y = (float)(gs.spacecraft[i].pos.y / SCALE);
        scaled_craft_pos[i].z = (float)(gs.spacecraft[i].pos.z / SCALE);
    }

    if (wp.draw_lines_between_bodies) {
        // draw lines between planets to show distance
        for (int i = 0; i < gb.count; i++) {
            const vec3_f pp1 = scaled_body_pos[i];
            int pp2idx = i + 1;
            if (i + 1 > gb.count - 1) pp2idx = 0;
            const vec3_f pp2 = scaled_body_pos[pp2idx];
            addLine(line_batch, pp1.x, pp1.y, pp1.z, pp2.x, pp2.y, pp2.z, 1, 1, 1);
        }
    }
    if (wp.draw_inclination_height) {
        for (int i = 0; i < gb.count; i++) {
            const vec3_f pp = scaled_body_pos[i];
            float r, g, b;
            if (pp.z > 0) { r = 0.5f, g = 0.5f; b = 1.0f; }
            else { r = 1.0f, g = 0.5f; b = 0.5f; }
            addLine(line_batch, pp.x, pp.y, pp.z, pp.x, pp.y, 0, r, g, b);
        }
    }

    // draw rotation axes for all bodies
    for (int i = 0; i < gb.count; i++) {
        const body_t body = gb.bodies[i];
        if (body.rotational_v != 0.0) {
            // get planet position (already scaled)
            const vec3_f planet_pos = scaled_body_pos[i];

            // extract rotation axis from attitude quaternion
            const vec3 z_axis = {0.0, 0.0, 1.0};
            const vec3 rotation_axis = quaternionRotate(body.attitude, z_axis);

            // scale the axis to extend beyond the planet
            const float axis_length = (float)(body.radius / SCALE) * 1.5f;
            const vec3_f axis_end1 = {
                planet_pos.x + (float)rotation_axis.x * axis_length,
                planet_pos.y + (float)rotation_axis.y * axis_length,
                planet_pos.z + (float)rotation_axis.z * axis_length
            };
            const vec3_f axis_end2 = {
                planet_pos.x - (float)rotation_axis.x * axis_length,
                planet_pos.y - (float)rotation_axis.y * axis_length,
                planet_pos.z - (float)rotation_axis.z * axis_length
            };

            // draw the rotation axis line
            addLine(line_batch, axis_end1.x, axis_end1.y, axis_end1.z,
                               axis_end2.x, axis_end2.y, axis_end2.z,
                               0.0f, 1.0f, 1.0f);
        }
    }

    //////////////////////////////////////////////////
    // craft orbital visuals
    //////////////////////////////////////////////////
    // render orbit property lines
    for (int i = 0; i < gs.count; i++) {
        spacecraft_t craft = gs.spacecraft[i];
        const vec3_f craft_pos = scaled_craft_pos[i];
        const vec3_f body_pos = scaled_body_pos[craft.SOI_planet_id];

        // line from planet to craft!
        addLine(line_batch, craft_pos.x, craft_pos.y, craft_pos.z, body_pos.x, body_pos.y, body_pos.z, 1.0f, 1.0f, 1.0f);
        // line to view angle from equatorial plane
        addLine(line_batch, craft_pos.x, craft_pos.y, body_pos.z, body_pos.x, body_pos.y, body_pos.z, 1.0f, 0.0f, 0.0f);
        // line to connect the triangle
        addLine(line_batch, craft_pos.x, craft_pos.y, craft_pos.z, craft_pos.x, craft_pos.y, body_pos.z, 0.0f, 1.0f, 0.0f);
    }

    renderCraftPaths(&sim, line_batch, craft_paths);
    renderPlanetPaths(&sim, line_batch, planet_paths);
}
