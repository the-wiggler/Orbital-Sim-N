#version 330 core
in vec2 uv;
out vec4 FragColor;
uniform sampler2D tex;
uniform vec3 color;

void main() {
    float a = texture(tex, uv).r;
    if (a < 0.5) discard;
    FragColor = vec4(color, 1.0);
}