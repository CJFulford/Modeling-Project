#include "Icon.h"
#include "lodepng.h"
#include "Shader.h"

using namespace std;
using namespace glm;


Icon::Icon(string file, vec2 pos)
{
	vertexArray = 0;
	vertexBuffer = 0;
	uvBuffer = 0;
	textureID = 0;
	
	filename = file;
	position = pos;

	imageHeight = 32;
	imageWidth = 32;
	program = generateProgram("Icon.vert", "Icon.frag");
	

	uvs.push_back(vec2(0.f, 1.f));
	uvs.push_back(vec2(1.f, 1.f));
	uvs.push_back(vec2(0.f, 0.f));
	uvs.push_back(vec2(1.f, 0.f));
}


void Icon::generateBuffer()
{
	glGenVertexArrays(1, &vertexArray);
	glBindVertexArray(vertexArray);

	glGenBuffers(1, &vertexBuffer);
	glGenBuffers(1, &uvBuffer);

}

void Icon::render()
{
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	glUseProgram(program);
	glBindVertexArray(vertexArray);
	texture.bind2DTexture(program, textureID, std::string("tex"));

	//draw 2 triangles 
	glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);

	glBindVertexArray(0);
	texture.unbind2DTexture();
	glUseProgram(0);

	glDisable(GL_BLEND);

}

void Icon::update()
{
	glBindVertexArray(vertexArray);

	verts.clear();	
	getTriangles();

	glBindBuffer(GL_ARRAY_BUFFER, vertexBuffer);
	glBufferData(GL_ARRAY_BUFFER, sizeof(vec2) * verts.size(), verts.data(), GL_STATIC_DRAW);
	glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 0, NULL);
	glEnableVertexAttribArray(0);

	glBindBuffer(GL_ARRAY_BUFFER, uvBuffer);
	glBufferData(GL_ARRAY_BUFFER, sizeof(vec2) * uvs.size(), uvs.data(), GL_STATIC_DRAW);
	glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 0, NULL);
	glEnableVertexAttribArray(1);

	glBindVertexArray(0);
}

//void Icon::update(string file)
//{
//	filename = file;
//}

void Icon::loadImages()
{
	getTriangles();
	generateBuffer();

	image.clear();
	glDeleteTextures(1, &textureID);

	unsigned int error = lodepng::decode(image, imageWidth, imageHeight, filename);

	if (error)
	{
		std::cout << "error " << error << ":" << lodepng_error_text(error) << std::endl;
	}
	//creating 2D texture
	textureID = texture.create2DTexture(image, imageWidth, imageHeight);

}

void Icon::getTriangles()
{
	verts.push_back(position + vec2(-.03f, -.03f));
	verts.push_back(position + vec2( .03f, -.03f));
	verts.push_back(position + vec2(-.03f,  .03f));
	verts.push_back(position + vec2( .03f,  .03f));
}