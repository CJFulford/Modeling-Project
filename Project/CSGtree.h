#include "Header.h"
#include "Icon.h"
#include "ShaderBuilder.h"
#include <vector>

#pragma once
class CSGtree
{
public:
	CSGtree();

    std::vector<Icon> icons;
	std::vector<glm::vec2> verts;
	std::vector<glm::vec3> colours;
    std::vector<glm::vec2> info; // <tree level, object type>



	void render();
    void constructInfo(Object *obj, int level);
    void update();

private:
	void generateBuffer();
    void makeVerts();
    GLuint vertexBuffer;
    GLuint colourBuffer;
    GLuint vertexArray;
    GLuint program;

};

