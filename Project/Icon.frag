#version 430 core

in vec2 UV;

uniform sampler2D tex;

out vec4 color;

void main(void)
{
    //color = vec4(1.f, 0.f, 0.f, 0.f); 
	vec4 precolor = texture(tex, UV);

	vec3 white = vec3(0.6,0.6,0.6);

	if(precolor.x >= white.x && precolor.y >= white.y && precolor.z >= white.z)
		precolor.w = 0.0;	

	color = precolor;
}