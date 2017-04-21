#include "Header.h"
#include "Icon.h"
#include <vector>

#pragma once
class CSGtree
{
public:
	CSGtree();
	void render();
    void constructInfo(Object * obj1, Object * obj2, int level);
    void constructInfo(Object *obj, int level);
    void update();

private:
    std::vector<Icon*> icons;
    std::vector<glm::vec2> verts;
    std::vector<glm::vec3> colours;
    std::vector<glm::vec2> info; // <tree level, object type>

    GLuint vertexBuffer;
    GLuint colourBuffer;
    GLuint vertexArray;
    GLuint program;

	void generateBuffer();
    void makeVerts();
};

