#include "Header.h"

#pragma once
class TList
{
public:
    TList();

    std::vector<glm::vec2> verts;

    GLuint vertexBuffer;
    GLuint vertexArray;
    GLuint program;

    void render();
    void getLines(Ray *ray);

private:
    void generateBuffer();
    void addToVerts(Ray *ray, int level);

};

