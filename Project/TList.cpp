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

    getLines(NULL);

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

    glLineWidth(5);
    glDrawArrays(GL_LINES, 10, verts.size());       //drawing thick lines over thin

    glBindVertexArray(0);
    glUseProgram(0);

}

void TList::addToVerts(Ray *ray, int level)
{
    float y = -0.1 + (level * -.2f);
    for (int i = 0; i < ray->volumes.size(); i++)
    {
        verts.push_back(vec2(ray->volumes[i].entrance, y));
        verts.push_back(vec2(ray->volumes[i].exit, y));
    }
    ray->volumes.clear();
}

void TList::getLines(Ray *ray)
{
    glBindVertexArray(vertexArray);

    verts.clear();
    float y = -0.1;
    for (int i = 0; i < 5; i++)
    {
        verts.push_back(vec2(-.9, y));
        verts.push_back(vec2(-.1, y));
        y -= 0.2;
    }


    // only one object selected
    if (selected1 != -1 && selected2 == -1)
    {
        Ray tempRay(ray->origin, ray->direction);
        objectVec[selected1]->getVolume(&tempRay);
    }
    // 2 objects selected
    else if (selected1 != -1 && selected2 != -1)
    {
        Ray tempRay(ray->origin, ray->direction);
        objectVec[selected1]->getVolume(&tempRay);
        addToVerts(&tempRay, 0);

        objectVec[selected2]->getVolume(&tempRay);
        addToVerts(&tempRay, 1);

        Union(objectVec[selected1], objectVec[selected2]).getVolume(&tempRay);
        addToVerts(&tempRay, 2);

        Intersection(objectVec[selected1], objectVec[selected2]).getVolume(&tempRay);
        addToVerts(&tempRay, 3);

        Difference(objectVec[selected1], objectVec[selected2]).getVolume(&tempRay);
        addToVerts(&tempRay, 4);

        objectVec[selected2]->differenceB = false;
    }



    glBindBuffer(GL_ARRAY_BUFFER, vertexBuffer);
    glBufferData(GL_ARRAY_BUFFER, sizeof(verts[0])*verts.size(), &verts[0], GL_STATIC_DRAW);
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 0, NULL);
    glEnableVertexAttribArray(0);

}