#version 430 core

layout (location = 0) in vec2 vertex;
layout (location = 1) in vec3 color;

out vec3 colour;

void main(void)
{
    colour = color;
    gl_Position = vec4(vertex, 0.f, 1.f);
}