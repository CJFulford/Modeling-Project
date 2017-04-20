#include "Header.h"
#include "texture.h"
#include <vector>
#include <string.h>

#pragma once
class Icon
{
public:
	Icon(std::string file, glm::vec2 pos);
	~Icon();

	glm::vec2 position;

	std::vector<glm::vec2> verts;
	std::vector<glm::vec2> uvs;

	GLuint vertexBuffer;
	GLuint vertexArray;
	GLuint uvBuffer;
	GLuint textureID;
	GLuint program;

	Texture texture;
	std::vector<unsigned char> image;
	unsigned int imageWidth;
	unsigned int imageHeight;

	std::string filename;

	void render();
	void update();

	void loadImages();

private:
	void generateBuffer();
	void getTriangles();

};

