#version 430 core

in vec2 UV;

uniform sampler2D tex;

out vec4 color;

void main(void)
{
    color = vec4(1.f, 0.f, 0.f, 0.f); 
	//vec4 color = texture(tex, UV);
}