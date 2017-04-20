#include "Header.h"
#include "Icon.h"
#include "ShaderBuilder.h"
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


	void render();
    void constructInfo(Object *obj, int level);

private:
	void generateBuffer();
    void update();
    void makeVerts();
    GLuint vertexBuffer;
    GLuint colourBuffer;
    GLuint vertexArray;
    GLuint program;

};

