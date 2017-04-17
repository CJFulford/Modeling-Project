#include "TList.h"
#include "Shader.h"

using namespace std;
using namespace glm;


TList::TList()
{
    vertexArray = 0;
    vertexBuffer = 0;
    program = generateProgram("general.vert", "general.frag");
    generateBuffer();
}


void TList::generateBuffer()
{
    glGenVertexArrays(1, &vertexArray);
    glBindVertexArray(vertexArray);

    getLines();

    glGenBuffers(1, &vertexBuffer);
    glBindBuffer(GL_ARRAY_BUFFER, vertexBuffer);
    glBufferData(GL_ARRAY_BUFFER, sizeof(verts[0]) * verts.size(), &verts[0], GL_STATIC_DRAW);
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 0, NULL);
    glEnableVertexAttribArray(0);

    glBindVertexArray(0);

}

void TList::render()
{
    glUseProgram(program);
    glBindVertexArray(vertexArray);

    glLineWidth(1);
    glDrawArrays(GL_LINES, 0, 10);                 //drawing five thin lines

    //glLineWidth(5);
   // glDrawArrays(GL_LINES, 10, verts.size());       //drawing thick lines over thin

    glBindVertexArray(0);
    glUseProgram(0);

}

void TList::getLines()
{
    glBindVertexArray(vertexArray);

    verts.clear();
    float y = -0.15;
    for (int i = 0; i < 5; i++)
    {
        verts.push_back(vec2(-.9, y));
        verts.push_back(vec2(-.1, y));
        y -= 0.125;
    }

    glBindBuffer(GL_ARRAY_BUFFER, vertexBuffer);
    glBufferData(GL_ARRAY_BUFFER, sizeof(verts[0])*verts.size(), &verts[0], GL_STATIC_DRAW);
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 0, NULL);
    glEnableVertexAttribArray(0);

}