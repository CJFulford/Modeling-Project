#include "Header.h"

#pragma once
class TList
{
public:
    TList();
    void render();
    void getLines(Ray *ray);

private:
    GLuint vertexBuffer;
    GLuint colourBuffer;
    GLuint vertexArray;
    GLuint program;
    void generateBuffer();
    void addToVerts(Ray *ray, int level, glm::vec3 colour);
    std::vector<glm::vec2> verts;
    std::vector<glm::vec3> colours;

};

