#include "Header.h"

#pragma once
class TList
{
public:
    TList();

    std::vector<glm::vec2> verts;
    std::vector<glm::vec3> colours;

    GLuint vertexBuffer;
    GLuint colourBuffer;
    GLuint vertexArray;
    GLuint program;

    void render();
    void getLines(Ray *ray);

private:
    void generateBuffer();
    void addToVerts(Ray *ray, int level, glm::vec3 colour);

};

