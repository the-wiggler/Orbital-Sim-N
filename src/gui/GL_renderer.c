//
// Created by java on 12/25/25.
//

#include "GL_renderer.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include <inttypes.h>
#include "models.h"

#include "../globals.h"
#include "../math/matrix.h"
#include "../types.h"

extern void displayError(const char* title, const char* message);

char* loadShaderSource(const char* filepath) {
    // use two alternating buffers to support loading vertex and fragment shaders simultaneously
    static char shader_buffer_0[SHADER_BUFFER_SIZE];
    static char shader_buffer_1[SHADER_BUFFER_SIZE];
    static int buffer_index = 0;

    // alternate between buffers
    char* shader_buffer = (buffer_index == 0) ? shader_buffer_0 : shader_buffer_1;
    buffer_index = 1 - buffer_index;

    #ifdef _WIN32
    FILE* file;
    fopen_s(&file, filepath, "r");
    #else
    FILE* file = fopen(filepath, "rb");
    #endif
    if (!file) {
        fprintf(stderr, "Failed to open shader file: %s\n", filepath);
        return NULL;
    }

    fseek(file, 0, SEEK_END);
    const long size = ftell(file);
    fseek(file, 0, SEEK_SET);

    // validate size
    if (size >= SHADER_BUFFER_SIZE) {
        fprintf(stderr, "Shader file too large: %ld bytes (max %d): %s\n", size, SHADER_BUFFER_SIZE, filepath);
        fclose(file);
        return NULL;
    }

    const size_t bytesRead = fread(shader_buffer, 1, size, file);
    shader_buffer[bytesRead] = '\0';  // Use actual bytes read, not file size

    fclose(file);
    return shader_buffer;
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
        return 0;
    }

    glDeleteShader(vertexShader);
    glDeleteShader(fragmentShader);

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

    // NOLINTNEXTLINE(performance-no-int-to-ptr) - required for OpenGL buffer offset
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)(uintptr_t)0);
    glEnableVertexAttribArray(0);

    // NOLINTNEXTLINE(performance-no-int-to-ptr) - required for OpenGL buffer offset
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)(uintptr_t)(3 * sizeof(float)));
    glEnableVertexAttribArray(1);

    // unbind
    glBindVertexArray(0);

    return vbo;
}

// initializes the openGL shader program with the .vert files we should have stored with the program. It also initializes
// relevant things used by the renderer
GLuint init_GL_shader() {
    // create the shader programs
    GLuint shaderProgram = createShaderProgram("shaders/simple.vert", "shaders/simple.frag");
    if (shaderProgram == 0) {
        displayError("Shader Error", "Failed to create shader program. Check console for details.");
        return 1;
    }

    // enable depth testing
    glEnable(GL_DEPTH_TEST);

    // enable blending for transparency
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    return shaderProgram;
}

// initializes all the assets we would want to use for the program, such as the VBOs and their models and the fonts.
GL_assets_t init_GL_assets(sim_properties_t* sim) {
    GL_assets_t assets = {0};

    // create buffer for cube shape
    assets.unit_cube = createVBO(UNIT_CUBE_VERTICES, sizeof(UNIT_CUBE_VERTICES));

    // create buffer for cone shape
    assets.cone = createVBO(CONE_VERTICES, sizeof(CONE_VERTICES));

    // create buffer for sphere shape
    sphere_mesh_t sphere_mesh;
    sphere_mesh = generateUnitSphere(15, 15);
    assets.sphere = createVBO(sphere_mesh.vertices, sphere_mesh.data_size);
    sim->window_params.planet_model_vertex_count = (int)sphere_mesh.vertex_count; // I couldn't think of a better way to do this ngl

    // create batch to hold all the line geometries we would ever want to draw
    assets.lines = createLineBatch(MAX_LINE_BATCH);

    // // craft path tracking
    // assets.craft_paths = {0};

    // initialize font for text rendering
    assets.font = initFont("assets/font.ttf", 24.0F);

    return assets;
}

void deleteVBO(const VBO_t vbo) {
    glDeleteVertexArrays(1, &vbo.VAO);
    glDeleteBuffers(1, &vbo.VBO);
}

// creates a view matrix that looks at target from cameraPos
mat4 createViewMatrix(const vec3_f cameraPos, const vec3_f target) {
    mat4 viewMatrix;
    vec3_f forward = vec3_f_normalize(vec3_f_sub(target, cameraPos));
    vec3_f up_vec = {0.0F, 0.0F, 1.0F}; // z is "UP"
    vec3_f right = vec3_f_normalize(vec3_f_cross(forward, up_vec));
    up_vec = vec3_f_cross(right, forward);

    // build view matrix
    viewMatrix.m[0] = right.x;
    viewMatrix.m[1] = up_vec.x;
    viewMatrix.m[2] = -forward.x;
    viewMatrix.m[3] = 0.0F;
    viewMatrix.m[4] = right.y;
    viewMatrix.m[5] = up_vec.y;
    viewMatrix.m[6] = -forward.y;
    viewMatrix.m[7] = 0.0F;
    viewMatrix.m[8] = right.z;
    viewMatrix.m[9] = up_vec.z;
    viewMatrix.m[10] = -forward.z;
    viewMatrix.m[11] = 0.0F;
    viewMatrix.m[12] = -vec3_f_dot(right, cameraPos);
    viewMatrix.m[13] = -vec3_f_dot(up_vec, cameraPos);
    viewMatrix.m[14] = vec3_f_dot(forward, cameraPos);
    viewMatrix.m[15] = 1.0F;
    return viewMatrix;
}

// creates a perspective projection matrix
// NOLINTNEXTLINE(bugprone-easily-swappable-parameters) - parameters are semantically different
mat4 createProjectionMatrix(const float fov, const float aspect, const float near_plane, const float far_plane) {
    mat4 projMatrix;
    const float focal_length = 1.0F / tanf(fov * 0.5F);

    projMatrix.m[0] = focal_length / aspect;
    projMatrix.m[1] = 0.0F;
    projMatrix.m[2] = 0.0F;
    projMatrix.m[3] = 0.0F;

    projMatrix.m[4] = 0.0F;
    projMatrix.m[5] = focal_length;
    projMatrix.m[6] = 0.0F;
    projMatrix.m[7] = 0.0F;

    projMatrix.m[8] = 0.0F;
    projMatrix.m[9] = 0.0F;
    projMatrix.m[10] = (far_plane + near_plane) / (near_plane - far_plane);
    projMatrix.m[11] = -1.0F;

    projMatrix.m[12] = 0.0F;
    projMatrix.m[13] = 0.0F;
    projMatrix.m[14] = (2.0F * far_plane * near_plane) / (near_plane - far_plane);
    projMatrix.m[15] = 0.0F;
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
        sim.window_params.cam_target.x + (sim.window_params.camera_pos.x * sim.window_params.zoom),
        sim.window_params.cam_target.y + (sim.window_params.camera_pos.y * sim.window_params.zoom),
        sim.window_params.cam_target.z + (sim.window_params.camera_pos.z * sim.window_params.zoom)
    };

    // create view matrix
    const mat4 viewMatrix = createViewMatrix(cameraPos, sim.window_params.cam_target);

    // create projection matrix
    const float aspect = sim.window_params.window_size_x / sim.window_params.window_size_y;
    const mat4 projMatrix = createProjectionMatrix(PI_OVER_4_f, aspect, 0.1F, 1000000000.0F);

    // set matrices in shader
    GL_setMatrixUniform(shaderProgram, "view", &viewMatrix);
    GL_setMatrixUniform(shaderProgram, "projection", &projMatrix);
}

// sphere mesh generation
sphere_mesh_t generateUnitSphere(const unsigned int stacks, const unsigned int sectors) {
    sphere_mesh_t sphere = {0};
    sphere.vertex_count = (size_t)stacks * sectors * 6;
    sphere.data_size = sphere.vertex_count * 6 * sizeof(float);

    // validate vertex count against static buffer size
    const size_t required_floats = sphere.vertex_count * 6;
    if (required_floats > MAX_SPHERE_VERTICES) {
        fprintf(stderr, "Sphere vertices exceed buffer: %zu floats required (max %d)\n", required_floats, MAX_SPHERE_VERTICES);
        sphere.vertex_count = 0;
        sphere.data_size = 0;
        return sphere;
    }

    float* data = sphere.vertices;

    // generate vertices for each stack and sector
    for (unsigned int i = 0; i < stacks; ++i) {
        const float theta1 = (float)i * PI_f / (float)stacks;
        const float theta2 = (float)(i + 1) * PI_f / (float)stacks;

        for (unsigned int j = 0; j < sectors; ++j) {
            const float phi1 = (float)j * 2.0F * PI_f / (float)sectors;
            const float phi2 = (float)(j + 1) * 2.0F * PI_f / (float)sectors;
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

// line rendering stuff
line_batch_t createLineBatch(const size_t max_lines) {
    line_batch_t batch = {0};

    // validate capacity against static buffer size
    if (max_lines > MAX_LINE_BATCH) {
        fprintf(stderr, "Line batch capacity %zu exceeds maximum %d\n", max_lines, MAX_LINE_BATCH);
        batch.capacity = MAX_LINE_BATCH;
    } else {
        batch.capacity = max_lines;
    }
    batch.count = 0;

    // calculate vertex buffer size
    const size_t vertex_array_size = batch.capacity * 2 * 6 * sizeof(float);

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
    // NOLINTNEXTLINE(performance-no-int-to-ptr) - required for OpenGL buffer offset
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)(uintptr_t)(3 * sizeof(float)));
    glEnableVertexAttribArray(1);

    glBindVertexArray(0);

    return batch;
}

// NOLINTNEXTLINE(bugprone-easily-swappable-parameters) - start/end points are semantically distinct
void addLine(line_batch_t* batch, const float start_x, const float start_y, const float start_z, const float end_x, const float end_y, const float end_z, const float red, const float green, const float blue) {
    if (!batch || batch->count >= batch->capacity) {
        return;
    }

    // calculate index for this line's data
    const size_t index = batch->count * 12; // 12 floats per line (2 vertices * 6 floats)

    // first vertex
    batch->vertices[index + 0] = start_x;
    batch->vertices[index + 1] = start_y;
    batch->vertices[index + 2] = start_z;
    batch->vertices[index + 3] = red;
    batch->vertices[index + 4] = green;
    batch->vertices[index + 5] = blue;

    // second vertex
    batch->vertices[index + 6] = end_x;
    batch->vertices[index + 7] = end_y;
    batch->vertices[index + 8] = end_z;
    batch->vertices[index + 9] = red;
    batch->vertices[index + 10] = green;
    batch->vertices[index + 11] = blue;

    batch->count++;
}

void renderLines(line_batch_t* batch, const GLuint shader_program) {
    if (!batch || batch->count == 0) {
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

#define STB_TRUETYPE_IMPLEMENTATION
#include "stb_truetype.h"

#define ATLAS 512
#define MAX_CHARS 4096

static stbtt_bakedchar cdata[96];

font_t initFont(const char* path, const float size) {
    font_t font_struct = {0};

    #ifdef _WIN32
    FILE* file_ptr;
    fopen_s(&file_ptr, path, "rb");
    #else
    FILE* file_ptr = fopen(path, "rb");
    #endif
    if (!file_ptr) {
        fprintf(stderr, "Failed to open font file: %s\n", path);
        return font_struct;
    }
    fseek(file_ptr, 0, SEEK_END);
    const long file_len = ftell(file_ptr);
    fseek(file_ptr, 0, SEEK_SET);
    unsigned char* font_data = malloc((size_t)file_len);
    fread(font_data, 1, (size_t)file_len, file_ptr);
    fclose(file_ptr);

    unsigned char* glyph_bitmap = malloc((size_t)ATLAS * ATLAS);
    // NOLINTNEXTLINE(clang-analyzer-security.ArrayBound) - false positive in stbtt_BakeFontBitmap
    stbtt_BakeFontBitmap(font_data, 0, size, glyph_bitmap, ATLAS, ATLAS, 32, 96, cdata);
    free(font_data);

    glGenTextures(1, &font_struct.tex);
    glBindTexture(GL_TEXTURE_2D, font_struct.tex);

    GLint internalFormat = GL_RED;

#ifdef __EMSCRIPTEN__
    // WebGL needs a different internal format
    internalFormat = GL_R8;
#endif

    glTexImage2D(GL_TEXTURE_2D, 0, internalFormat, ATLAS, ATLAS, 0, GL_RED, GL_UNSIGNED_BYTE, glyph_bitmap);

    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    free(glyph_bitmap);

    font_struct.shader = createShaderProgram("shaders/text.vert", "shaders/text.frag");
    if (font_struct.shader == 0) {
        fprintf(stderr, "Failed to create text shader program\n");
        glDeleteTextures(1, &font_struct.tex);
        return font_struct;
    }

    glGenVertexArrays(1, &font_struct.vao);
    glGenBuffers(1, &font_struct.vbo);
    glBindVertexArray(font_struct.vao);
    glBindBuffer(GL_ARRAY_BUFFER, font_struct.vbo);
    glBufferData(GL_ARRAY_BUFFER, (long long)MAX_CHARS * 24 * sizeof(float), NULL, GL_DYNAMIC_DRAW);
    glVertexAttribPointer(0, 4, GL_FLOAT, GL_FALSE, 4 * sizeof(float), 0);
    glEnableVertexAttribArray(0);

    printf("Successfully initialized font: %s\n", path);
    return font_struct;
}

void addText(font_t* font, float pos_x, const float pos_y, const char* text, const float scale) {
    for (; *text; text++) {
        const int char_index = *text - 32;
        if (char_index < 0 || char_index >= 96) {
            continue;
        }
        const stbtt_bakedchar* glyph_char = &cdata[char_index];

        float glyph_x0 = pos_x + (glyph_char->xoff * scale);
        float glyph_y0 = pos_y + (glyph_char->yoff * scale);
        float glyph_x1 = glyph_x0 + ((float)(glyph_char->x1 - glyph_char->x0) * scale);
        float glyph_y1 = glyph_y0 + ((float)(glyph_char->y1 - glyph_char->y0) * scale);
        float tex_s0 = (float)glyph_char->x0 / (float)ATLAS;
        float tex_t0 = (float)glyph_char->y0 / (float)ATLAS;
        float tex_s1 = (float)glyph_char->x1 / (float)ATLAS;
        float tex_t1 = (float)glyph_char->y1 / (float)ATLAS;

        float* vertex_data = &font->verts[(size_t)font->count++ * 24];
        vertex_data[0]=glyph_x0; vertex_data[1]=glyph_y0; vertex_data[2]=tex_s0; vertex_data[3]=tex_t0;
        vertex_data[4]=glyph_x1; vertex_data[5]=glyph_y0; vertex_data[6]=tex_s1; vertex_data[7]=tex_t0;
        vertex_data[8]=glyph_x1; vertex_data[9]=glyph_y1; vertex_data[10]=tex_s1; vertex_data[11]=tex_t1;
        vertex_data[12]=glyph_x0; vertex_data[13]=glyph_y0; vertex_data[14]=tex_s0; vertex_data[15]=tex_t0;
        vertex_data[16]=glyph_x1; vertex_data[17]=glyph_y1; vertex_data[18]=tex_s1; vertex_data[19]=tex_t1;
        vertex_data[20]=glyph_x0; vertex_data[21]=glyph_y1; vertex_data[22]=tex_s0; vertex_data[23]=tex_t1;

        pos_x += glyph_char->xadvance * scale;
    }
}

// NOLINTNEXTLINE(bugprone-easily-swappable-parameters) - window dimensions are semantically distinct from colors
void renderText(font_t* font, const float window_width, const float window_height, const float red, const float green, const float blue) {
    if (!font->count) {
        return;
    }

    glDisable(GL_DEPTH_TEST);
    glUseProgram(font->shader);

    const float proj[16] = {
        2.0F/window_width, 0.0F, 0.0F, 0.0F,
        0.0F, -2.0F/window_height, 0.0F, 0.0F,
        0.0F, 0.0F, -1.0F, 0.0F,
        -1.0F, 1.0F, 0.0F, 1.0F
    };
    glUniformMatrix4fv(glGetUniformLocation(font->shader, "proj"), 1, GL_FALSE, proj);
    glUniform3f(glGetUniformLocation(font->shader, "color"), red, green, blue);

    glBindTexture(GL_TEXTURE_2D, font->tex);
    glBindBuffer(GL_ARRAY_BUFFER, font->vbo);
    glBufferSubData(GL_ARRAY_BUFFER, 0, (GLsizeiptr)((size_t)font->count * 24 * sizeof(float)), font->verts);
    glBindVertexArray(font->vao);
    glDrawArrays(GL_TRIANGLES, 0, font->count * 6);

    glEnable(GL_DEPTH_TEST);
    font->count = 0;
}








// render the coordinate plane to the screen
void renderCoordinatePlane(const sim_properties_t sim, line_batch_t* line_batch) {
    const float scale = sim.window_params.zoom;

    // X axis (red) - horizontal
    addLine(line_batch, -10.0F * scale, 0.0F, 0.0F, 10.0F * scale, 0.0F, 0.0F, 0.3F, 0.0F, 0.0F);

    // Y axis (green) - horizontal
    addLine(line_batch, 0.0F, -10.0F * scale, 0.0F, 0.0F, 10.0F * scale, 0.0F, 0.0F, 0.3F, 0.0F);

    // Z axis (blue) - vertical "up"
    addLine(line_batch, 0.0F, 0.0F, -10.0F * scale, 0.0F, 0.0F, 10.0F * scale, 0.0F, 0.0F, 0.3F);

    // perspective lines in X-Y plane (gray)
    addLine(line_batch, 10.0F * scale, 10.0F * scale, 0.0F, -10.0F * scale, -10.0F * scale, 0.0F, 0.3F, 0.3F, 0.3F);
    addLine(line_batch, 10.0F * scale, -10.0F * scale, 0.0F, -10.0F * scale, 10.0F * scale, 0.0F, 0.3F, 0.3F, 0.3F);
}

// render the sim planets to the screen
void renderPlanets(const sim_properties_t sim, const GLuint shader_program, const VBO_t planet_shape_buffer) {
    glBindVertexArray(planet_shape_buffer.VAO);
    for (int i = 0; i < sim.global_bodies.count; i++) {
        const body_t* body = &sim.global_bodies.bodies[i];

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
        glDrawArrays(GL_TRIANGLES, 0, sim.window_params.planet_model_vertex_count);
    }

    // draw sphere of influence
    if (sim.window_params.draw_planet_SOI) {
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        glDepthMask(GL_FALSE);
        glUniform1i(glGetUniformLocation(shader_program, "useOverride"), 1);
        glUniform4f(glGetUniformLocation(shader_program, "colorOverride"), 0.0F, 0.5F, 1.0F, 0.2F);
        glBindVertexArray(planet_shape_buffer.VAO);
        for (int i = 0; i < sim.global_bodies.count; i++) {
            const body_t* body = &sim.global_bodies.bodies[i];
            const float size_scale_factor = (float)body->SOI_radius / SCALE;
            const mat4 scale_mat = mat4_scale(size_scale_factor, size_scale_factor, size_scale_factor);
            const mat4 rotation_mat = quaternionToMatrix(body->attitude);
            const mat4 translate_mat = mat4_translation((float)body->pos.x / SCALE,(float)body->pos.y / SCALE,(float)body->pos.z / SCALE);
            const mat4 temp = mat4_mul(rotation_mat, scale_mat);
            mat4 planet_model = mat4_mul(translate_mat, temp);
            GL_setMatrixUniform(shader_program, "model", &planet_model);
            glDrawArrays(GL_TRIANGLES, 0, sim.window_params.planet_model_vertex_count);
        }
        glUniform1i(glGetUniformLocation(shader_program, "useOverride"), 0);
        glDepthMask(GL_TRUE);
        glDisable(GL_BLEND);
    }
}

// render the sim crafts to the renderer
void renderCrafts(const sim_properties_t sim, const GLuint shader_program, const VBO_t craft_shape_buffer) {
    glBindVertexArray(craft_shape_buffer.VAO);
    for (int i = 0; i < sim.global_spacecraft.count; i++) {
        const spacecraft_t* craft = &sim.global_spacecraft.spacecraft[i];

        // create a scale matrix
        const float size_scale_factor = 0.05F;
        const mat4 scale_matrix = mat4_scale(size_scale_factor, size_scale_factor * 2, size_scale_factor);
        const mat4 rotation_matrix = quaternionToMatrix(craft->attitude);
        const mat4 translation_matrix = mat4_translation(
            (float)craft->pos.x / SCALE,
            (float)craft->pos.y / SCALE,
            (float)craft->pos.z / SCALE);

        const mat4 temp = mat4_mul(rotation_matrix, scale_matrix);
        mat4 model = mat4_mul(translation_matrix, temp);

        // apply matrix and render to screen
        GL_setMatrixUniform(shader_program, "model", &model);
        glDrawArrays(GL_TRIANGLES, 0, 48);
    }
}

// render the stats on the screen
void renderStats(const sim_properties_t sim, font_t* font) {

    // calculate proper line height
    const float line_height = 20.0F;
    const float cursor_starting_pos[2] = { 10.0F, line_height + 10.0F};
    float cursor_pos[2] = { cursor_starting_pos[0], cursor_starting_pos[1] };

    // text buffer used to hold text written to the stats window
    char text_buffer[64];

    // paused indication
    if (sim.window_params.sim_running) {
        snprintf(text_buffer, sizeof(text_buffer), "Sim running");
    } else {
        snprintf(text_buffer, sizeof(text_buffer), "Sim paused");
    }
    addText(font, cursor_pos[0], cursor_pos[1], text_buffer, 0.8F);
    cursor_pos[1] += line_height;

    // write time step
    snprintf(text_buffer, sizeof(text_buffer), "Step: %.4g", sim.window_params.time_step);
    addText(font, cursor_pos[0], cursor_pos[1], text_buffer, 0.8F);
    cursor_pos[1] += line_height;

    // time indication
    const double time = sim.window_params.sim_time / 3600;
    if (time < 72.0) {
        snprintf(text_buffer, sizeof(text_buffer), "Time: %.2F hrs", time);
    } else if (time < 8766.0) {
        snprintf(text_buffer, sizeof(text_buffer), "Time: %.2F days", time / 24);
    }
    addText(font, cursor_pos[0], cursor_pos[1], text_buffer, 0.8F);
    cursor_pos[1] += line_height;

    // spacer
    cursor_pos[1] += line_height;

    for (int i = 0; i < sim.global_spacecraft.count; i++) {
        const int closest_id = sim.global_spacecraft.spacecraft[i].orbital_elements.closest_planet_id;
        const int soi_id = sim.global_spacecraft.spacecraft[i].orbital_elements.SOI_planet_id;

        snprintf(text_buffer, sizeof(text_buffer), "%s", sim.global_spacecraft.spacecraft[i].name);
        addText(font, cursor_pos[0], cursor_pos[1], text_buffer, 0.8F);
        cursor_pos[1] += line_height;

        snprintf(text_buffer, sizeof(text_buffer), "Closest Planet: %s", sim.global_bodies.bodies[closest_id].name);
        addText(font, cursor_pos[0], cursor_pos[1], text_buffer, 0.7F);
        cursor_pos[1] += line_height;

        snprintf(text_buffer, sizeof(text_buffer), "Distance: %.2F km", (sqrt(sim.global_spacecraft.spacecraft[i].orbital_elements.closest_r_squared) / 1000.0) - (sim.global_bodies.bodies[closest_id].radius / 1000.0));
        addText(font, cursor_pos[0], cursor_pos[1], text_buffer, 0.7F);
        cursor_pos[1] += line_height;

        snprintf(text_buffer, sizeof(text_buffer), "In SOI of: %s", sim.global_bodies.bodies[soi_id].name);
        addText(font, cursor_pos[0], cursor_pos[1], text_buffer, 0.7F);
        cursor_pos[1] += line_height;

        snprintf(text_buffer, sizeof(text_buffer), "Semi Major Axis: %.4g m", sim.global_spacecraft.spacecraft[i].orbital_elements.semi_major_axis);
        addText(font, cursor_pos[0], cursor_pos[1], text_buffer, 0.7F);
        cursor_pos[1] += line_height;
        snprintf(text_buffer, sizeof(text_buffer), "Eccentricity: %.6F", sim.global_spacecraft.spacecraft[i].orbital_elements.eccentricity);
        addText(font, cursor_pos[0], cursor_pos[1], text_buffer, 0.7F);
        cursor_pos[1] += line_height;
        snprintf(text_buffer, sizeof(text_buffer), "Inclination: %.4F rad", sim.global_spacecraft.spacecraft[i].orbital_elements.inclination);
        addText(font, cursor_pos[0], cursor_pos[1], text_buffer, 0.7F);
        cursor_pos[1] += line_height;
        snprintf(text_buffer, sizeof(text_buffer), "Ascending Node: %.4F rad", sim.global_spacecraft.spacecraft[i].orbital_elements.ascending_node);
        addText(font, cursor_pos[0], cursor_pos[1], text_buffer, 0.7F);
        cursor_pos[1] += line_height;
        snprintf(text_buffer, sizeof(text_buffer), "Arg of Periapsis: %.4F rad", sim.global_spacecraft.spacecraft[i].orbital_elements.arg_periapsis);
        addText(font, cursor_pos[0], cursor_pos[1], text_buffer, 0.7F);
        cursor_pos[1] += line_height;
        snprintf(text_buffer, sizeof(text_buffer), "True Anomaly: %.4F rad", sim.global_spacecraft.spacecraft[i].orbital_elements.true_anomaly);
        addText(font, cursor_pos[0], cursor_pos[1], text_buffer, 0.7F);
        cursor_pos[1] += line_height;
        snprintf(text_buffer, sizeof(text_buffer), "Specific Energy: %.4F J", sim.global_spacecraft.spacecraft[i].orbital_elements.specific_E);
        addText(font, cursor_pos[0], cursor_pos[1], text_buffer, 0.7F);
        cursor_pos[1] += line_height * 1.5F;

        snprintf(text_buffer, sizeof(text_buffer), "Fuel: %.1F kg", sim.global_spacecraft.spacecraft[i].fuel_mass);
        addText(font, cursor_pos[0], cursor_pos[1], text_buffer, 0.7F);
        cursor_pos[1] += line_height;
        snprintf(text_buffer, sizeof(text_buffer), "Engine: %s", (sim.global_spacecraft.spacecraft[i].engine_on) ? "on" : "off");
        addText(font, cursor_pos[0], cursor_pos[1], text_buffer, 0.7F);
        cursor_pos[1] += line_height;
        snprintf(text_buffer, sizeof(text_buffer), "Absolute v: %f", vec3_mag(sim.global_spacecraft.spacecraft[i].vel));
        addText(font, cursor_pos[0], cursor_pos[1], text_buffer, 0.7F);
        cursor_pos[1] += line_height * 1.5F;

        snprintf(text_buffer, sizeof(text_buffer), "Target Planet: %s", sim.global_bodies.bodies[sim.global_spacecraft.spacecraft[i].target_body_id].name);
        addText(font, cursor_pos[0], cursor_pos[1], text_buffer, 0.7F);

        cursor_pos[1] += line_height;
    }
}

void renderCraftPaths(const sim_properties_t* sim, line_batch_t* line_batch, craft_path_storage_t* craft_paths) {
    window_params_t window_params = sim->window_params;
    spacecraft_properties_t global_spacecraft = sim->global_spacecraft;

    if (window_params.draw_craft_path && craft_paths->num_objects > 0) {
        // draw orbital paths for all crafts
        for (int craft_idx = 0; craft_idx < craft_paths->num_objects; craft_idx++) {
            const int base = craft_idx * craft_paths->capacity;
            for (int i = 1; i < craft_paths->counts[craft_idx]; i++) {
                const vec3 prev_pos = craft_paths->positions[base + i - 1];
                const vec3 curr_pos = craft_paths->positions[base + i];
                addLine(line_batch, (float)prev_pos.x, (float)prev_pos.y, (float)prev_pos.z,
                                    (float)curr_pos.x, (float)curr_pos.y, (float)curr_pos.z, 1.0F, 1.0F, 0.5F);
            }
        }
    }

    // initialize or resize craft paths if needed
    if (global_spacecraft.count > 0 && craft_paths->num_objects != global_spacecraft.count) {
        craft_paths->num_objects = global_spacecraft.count;
        craft_paths->capacity = PATH_CAPACITY;
        // clear counts for all spacecraft
        for (int i = 0; i < craft_paths->num_objects && i < MAX_SPACECRAFT; i++) {
            craft_paths->counts[i] = 0;
        }
    }

    // record craft paths
    if (global_spacecraft.count > 0) {
        for (int craft_idx = 0; craft_idx < global_spacecraft.count; craft_idx++) {
            const spacecraft_t* craft = &sim->global_spacecraft.spacecraft[craft_idx];
            const int idx = (craft_idx * craft_paths->capacity) + craft_paths->counts[craft_idx];
            if (craft_paths->counts[craft_idx] < craft_paths->capacity) {
                craft_paths->positions[idx].x = craft->pos.x / SCALE;
                craft_paths->positions[idx].y = craft->pos.y / SCALE;
                craft_paths->positions[idx].z = craft->pos.z / SCALE;
                craft_paths->counts[craft_idx]++;
            } else {
                const int base = craft_idx * craft_paths->capacity;
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

void renderPredictedOrbits(sim_properties_t sim, line_batch_t* line_batch) {
    // render predicted orbits for bodies (skip central body at index 0)
    for (int i = 1; i < sim.global_bodies.count; i++) {
        body_t body = sim.global_bodies.bodies[i];
        if (body.oe.SOI_planet_id == i) {
            continue;
        }

        vec3_f parent_pos = {
            (float)(sim.global_bodies.bodies[body.oe.SOI_planet_id].pos.x / SCALE),
            (float)(sim.global_bodies.bodies[body.oe.SOI_planet_id].pos.y / SCALE),
            (float)(sim.global_bodies.bodies[body.oe.SOI_planet_id].pos.z / SCALE)
        };

        float semi_latus_rectum = (float)(body.oe.semi_major_axis * (1 - body.oe.eccentricity * body.oe.eccentricity));

        mat4 rotation_z_0 = mat4_rotationZ(-(float)body.oe.ascending_node);
        mat4 rotation_x_i = mat4_rotationX(-(float)body.oe.inclination);
        mat4 rotation_z_w = mat4_rotationZ(-(float)body.oe.arg_periapsis);
        mat4 rotation_matrix = mat4_mul(rotation_z_0, mat4_mul(rotation_x_i, rotation_z_w));

        if (body.oe.eccentricity >= 1.0) {
            continue;
        }

        int path_res = 100;
        float periapsis_distance = semi_latus_rectum / (1 + (float)body.oe.eccentricity);
        vec3_f position_perifocal_0 = { periapsis_distance, 0, 0 };
        vec3_f prev_position_eci = vec3_transformByMat4(rotation_matrix, position_perifocal_0);

        for (int j = 1; j <= path_res; j++) {
            // adaptive sampling: sample eccentric anomaly uniformly for better distribution on high-e orbits
            float eccentric_anomaly = TWO_PI_f * ((float)j / (float)path_res);

            // convert eccentric anomaly to true anomaly
            float true_anomaly = 2.0f * atanf(
                sqrtf((1.0f + (float)body.oe.eccentricity) /
                      (1.0f - (float)body.oe.eccentricity)) *
                tanf(eccentric_anomaly / 2.0f)
            );

            float orbital_radius = semi_latus_rectum / (1 + (float)body.oe.eccentricity * cosf(true_anomaly));
            vec3_f position_perifocal = {
                orbital_radius * cosf(true_anomaly),
                orbital_radius * sinf(true_anomaly),
                0
            };
            vec3_f position_eci = vec3_transformByMat4(rotation_matrix, position_perifocal);

            addLine(line_batch,
                parent_pos.x + (prev_position_eci.x / SCALE), parent_pos.y + (prev_position_eci.y / SCALE), parent_pos.z + (prev_position_eci.z / SCALE),
                parent_pos.x + (position_eci.x / SCALE), parent_pos.y + (position_eci.y / SCALE), parent_pos.z + (position_eci.z / SCALE),
                0.5F, 0.5F, 0.5F);

            prev_position_eci = position_eci;
        }
    }

    // render predicted orbits for spacecraft
    if (sim.global_spacecraft.count > 0) {
        for (int i = 0; i < sim.global_spacecraft.count; i++) {
            spacecraft_t craft = sim.global_spacecraft.spacecraft[i];
            vec3_f body_pos = {
                (float)(sim.global_bodies.bodies[craft.orbital_elements.SOI_planet_id].pos.x / SCALE),
                (float)(sim.global_bodies.bodies[craft.orbital_elements.SOI_planet_id].pos.y / SCALE),
                (float)(sim.global_bodies.bodies[craft.orbital_elements.SOI_planet_id].pos.z / SCALE)
            };

            float semi_latus_rectum = (float)(craft.orbital_elements.semi_major_axis * (1 - craft.orbital_elements.eccentricity * craft.orbital_elements.eccentricity));

            mat4 rotation_z_0 = mat4_rotationZ(-(float)craft.orbital_elements.ascending_node);
            mat4 rotation_x_i = mat4_rotationX(-(float)craft.orbital_elements.inclination);
            mat4 rotation_z_w = mat4_rotationZ(-(float)craft.orbital_elements.arg_periapsis);
            mat4 rotation_matrix = mat4_mul(rotation_z_0, mat4_mul(rotation_x_i, rotation_z_w));

            if (craft.orbital_elements.eccentricity >= 1.0) {
                continue;
            }

            // increment through eccentric anomaly values to build the orbit with given parameters
            // adaptive sampling: eccentric anomaly provides better distribution on high-e orbits
            int path_res = 100;
            float periapsis_distance = semi_latus_rectum / (1 + (float)craft.orbital_elements.eccentricity);
            vec3_f position_perifocal_0 = { periapsis_distance, 0, 0 };
            vec3_f prev_position_eci = vec3_transformByMat4(rotation_matrix, position_perifocal_0);

            for (int j = 1; j <= path_res; j++) {
                // sample eccentric anomaly uniformly for better distribution on high-e orbits
                float eccentric_anomaly = TWO_PI_f * ((float)j / (float)path_res);

                // convert eccentric anomaly to true anomaly
                float true_anomaly = 2.0f * atanf(
                    sqrtf((1.0f + (float)craft.orbital_elements.eccentricity) /
                          (1.0f - (float)craft.orbital_elements.eccentricity)) *
                    tanf(eccentric_anomaly / 2.0f)
                );

                float orbital_radius = semi_latus_rectum / (1 + (float)craft.orbital_elements.eccentricity * cosf(true_anomaly));
                vec3_f position_perifocal = {
                    orbital_radius * cosf(true_anomaly),
                    orbital_radius * sinf(true_anomaly),
                    0
                };
                vec3_f position_eci = vec3_transformByMat4(rotation_matrix, position_perifocal);

                addLine(line_batch,
                    body_pos.x + (prev_position_eci.x / SCALE), body_pos.y + (prev_position_eci.y / SCALE), body_pos.z + (prev_position_eci.z / SCALE),
                    body_pos.x + (position_eci.x / SCALE), body_pos.y + (position_eci.y / SCALE), body_pos.z + (position_eci.z / SCALE),
                    0.0F, 0.0F, 1.0F);

                prev_position_eci = position_eci;
            }

        }
    }
}

// renders debug features when they are enabled
void renderVisuals(sim_properties_t sim, line_batch_t* line_batch, craft_path_storage_t* craft_paths) {
    window_params_t window_params = sim.window_params;
    spacecraft_properties_t global_spacecraft = sim.global_spacecraft;
    body_properties_t global_bodies = sim.global_bodies;

    // create temporary arrays for scaled positions
    vec3_f scaled_body_pos[global_bodies.count];
    vec3_f scaled_craft_pos[global_spacecraft.count];

    for (int i = 0; i < global_bodies.count; i++) {
        scaled_body_pos[i].x = (float)(global_bodies.bodies[i].pos.x / SCALE);
        scaled_body_pos[i].y = (float)(global_bodies.bodies[i].pos.y / SCALE);
        scaled_body_pos[i].z = (float)(global_bodies.bodies[i].pos.z / SCALE);
    }
    for (int i = 0; i < global_spacecraft.count; i++) {
        scaled_craft_pos[i].x = (float)(global_spacecraft.spacecraft[i].pos.x / SCALE);
        scaled_craft_pos[i].y = (float)(global_spacecraft.spacecraft[i].pos.y / SCALE);
        scaled_craft_pos[i].z = (float)(global_spacecraft.spacecraft[i].pos.z / SCALE);
    }

    if (window_params.draw_lines_between_bodies) {
        // draw lines between planets to show distance
        for (int i = 0; i < global_bodies.count; i++) {
            const vec3_f body_pos_1 = scaled_body_pos[i];
            int other_body_index = i + 1;
            if (i + 1 > global_bodies.count - 1) {
                other_body_index = 0;
            }
            const vec3_f body_pos_2 = scaled_body_pos[other_body_index];
            addLine(line_batch, body_pos_1.x, body_pos_1.y, body_pos_1.z, body_pos_2.x, body_pos_2.y, body_pos_2.z, 1.0F, 1.0F, 1.0F);
        }
    }
    if (window_params.draw_inclination_height) {
        for (int i = 0; i < global_bodies.count; i++) {
            const vec3_f body_scaled_pos = scaled_body_pos[i];
            float red = 1.0F;
            float green = 0.5F;
            float blue = 0.5F;
            if (body_scaled_pos.z > 0) {
                red = 0.5F;
                green = 0.5F;
                blue = 1.0F;
            }
            addLine(line_batch, body_scaled_pos.x, body_scaled_pos.y, body_scaled_pos.z, body_scaled_pos.x, body_scaled_pos.y, 0.0F, red, green, blue);
        }
    }

    // draw rotation axes for all bodies
    for (int i = 0; i < global_bodies.count; i++) {
        const body_t body = global_bodies.bodies[i];
        if (body.rotational_v != 0.0) {
            // get planet position (already scaled)
            const vec3_f planet_pos = scaled_body_pos[i];

            // extract rotation axis from attitude quaternion
            const vec3 z_axis = {0.0, 0.0, 1.0};
            const vec3 rotation_axis = quaternionRotate(body.attitude, z_axis);

            // scale the axis to extend beyond the planet
            const float axis_length = (float)(body.radius / SCALE) * 1.5F;
            const vec3_f axis_end1 = {
                planet_pos.x + ((float)rotation_axis.x * axis_length),
                planet_pos.y + ((float)rotation_axis.y * axis_length),
                planet_pos.z + ((float)rotation_axis.z * axis_length)
            };
            const vec3_f axis_end2 = {
                planet_pos.x - ((float)rotation_axis.x * axis_length),
                planet_pos.y - ((float)rotation_axis.y * axis_length),
                planet_pos.z - ((float)rotation_axis.z * axis_length)
            };

            // draw the rotation axis line
            addLine(line_batch, axis_end1.x, axis_end1.y, axis_end1.z,
                               axis_end2.x, axis_end2.y, axis_end2.z,
                               0.0F, 1.0F, 1.0F);
        }
    }

    //////////////////////////////////////////////////
    // craft orbital visuals
    //////////////////////////////////////////////////
    // render orbit property lines
    for (int i = 0; i < global_spacecraft.count; i++) {
        spacecraft_t craft = global_spacecraft.spacecraft[i];
        const vec3_f craft_pos = scaled_craft_pos[i];
        const vec3_f body_pos = scaled_body_pos[craft.orbital_elements.SOI_planet_id];
        const body_t planet = global_bodies.bodies[craft.orbital_elements.SOI_planet_id];

        // get planet's rotation axis (equatorial plane normal) and project spacecraft onto equatorial plane
        const vec3 rotation_axis_norm = vec3_normalize(quaternionRotate(planet.attitude, (vec3){0.0, 0.0, 1.0}));
        const vec3 rel_pos = vec3_sub(craft.pos, planet.pos);
        const vec3 projected_rel = vec3_sub(rel_pos, vec3_scale(rotation_axis_norm, vec3_dot(rel_pos, rotation_axis_norm)));
        const vec3_f projected_scaled = {
            (float)((planet.pos.x + projected_rel.x) / SCALE),
            (float)((planet.pos.y + projected_rel.y) / SCALE),
            (float)((planet.pos.z + projected_rel.z) / SCALE)
        };

        // line from planet to craft
        addLine(line_batch, craft_pos.x, craft_pos.y, craft_pos.z, body_pos.x, body_pos.y, body_pos.z, 1.0F, 1.0F, 1.0F);
        // line from projected point on equatorial plane to planet center
        addLine(line_batch, projected_scaled.x, projected_scaled.y, projected_scaled.z, body_pos.x, body_pos.y, body_pos.z, 1.0F, 0.0F, 0.0F);
        // line from craft to projected point (inclination height line)
        addLine(line_batch, craft_pos.x, craft_pos.y, craft_pos.z, projected_scaled.x, projected_scaled.y, projected_scaled.z, 0.0F, 1.0F, 0.0F);

        // draw line to auto burn target position
        for (int j = 0; j < craft.num_burns; j++) {
            if (craft.burn_properties[j].auto_burn) {
                const vec3_f target_pos = {
                    (float)(craft.burn_properties[j].auto_burn_final_pos.x / SCALE),
                    (float)(craft.burn_properties[j].auto_burn_final_pos.y / SCALE),
                    (float)(craft.burn_properties[j].auto_burn_final_pos.z / SCALE)
                };
                addLine(line_batch, craft_pos.x, craft_pos.y, craft_pos.z, target_pos.x, target_pos.y, target_pos.z, 0.0F, 1.0F, 1.0F);
            }
        }
    }

    renderCraftPaths(&sim, line_batch, craft_paths);

    renderPredictedOrbits(sim, line_batch);
}
