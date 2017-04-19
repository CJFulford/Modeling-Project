#version 430 core

layout (location = 0) in vec2 vertex;
layout (location = 1) in vec2 uvs;

out vec2 UV;

void main(void)
{
	UV = uvs;
    gl_Position = vec4(vertex, 0.f, 1.f);
}