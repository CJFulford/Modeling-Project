#include "Header.h"
#include "Icon.h"
#include "Shader.h"
#include <vector>

#pragma once
class CSGtree
{
public:
	CSGtree();

	std::vector<glm::vec2> verts;
	std::vector<glm::vec3> colours, tempColours;
	std::vector<glm::vec2> info;			//index = level;
									//int at index = # nodes on level
	std::vector<Icon> icons;
	std::vector<Icon> tempIcons;

	GLuint vertexBuffer;
	GLuint colourBuffer;
	GLuint vertexArray;
	GLuint program;

	void render();
	void update();
	void constructInfo(Object *obj, int level);
	void makeVerts();

private:
	void generateBuffer();

};

