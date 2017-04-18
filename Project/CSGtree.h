#include "Header.h"
#include <vector>

#pragma once
class CSGtree
{
public:
	CSGtree();

	std::vector<glm::vec2> verts;
	std::vector<int> info;			//index = level;
									//int at index = # nodes on level
	GLuint vertexBuffer;
	GLuint colourBuffer;
	GLuint vertexArray;
	GLuint program;

	void render();
	void update();
	void constructInfo(Object *obj);

private:
	void generateBuffer();

};
