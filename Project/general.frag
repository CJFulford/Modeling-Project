#version 430 core

out vec4 color;

in vec3 colour;

void main(void)
{
    color = vec4(colour, 0.f);
}