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

    void generateBuffer();
    void render();
    void getLines();

private:

};

