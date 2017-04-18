#include "TList.h"
#include "Shader.h"

using namespace std;
using namespace glm;


TList::TList()
{
    vertexArray = 0;
    vertexBuffer = 0;
    program = generateProgram("TList.vert", "TList.frag");
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

    glGenBuffers(1, &colourBuffer);
    glBindBuffer(GL_ARRAY_BUFFER, colourBuffer);
    glBufferData(GL_ARRAY_BUFFER, sizeof(colours[0]) * colours.size(), &colours[0], GL_STATIC_DRAW);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, NULL);
    glEnableVertexAttribArray(1);
    glBindVertexArray(0);
}

void TList::render()
{
    glUseProgram(program);
    glBindVertexArray(vertexArray);

    glLineWidth(3);
    glDrawArrays(GL_LINES, 0, 10);                 //drawing five thin lines

    if (verts.size() > 10)
    {
        glLineWidth(10);
        glDrawArrays(GL_LINES, 10, verts.size() - 10);       //drawing thick lines over thin
    }

    glBindVertexArray(0);
    glUseProgram(0);

}

void TList::addToVerts(Ray *ray, int level, vec3 colour)
{
    float y = -0.5 + (level * -.1f);
    for (int i = 0; i < ray->volumes.size(); i++)
    {
        verts.push_back(vec2(ray->volumes[i].entrance, y));
        verts.push_back(vec2(ray->volumes[i].exit, y));
        colours.push_back(colour);
        colours.push_back(colour);
    }
    ray->volumes.clear();
}

void TList::getLines(Ray *ray)
{
    verts.clear();
    colours.clear();
    float y = -0.5;
    for (int i = 0; i < 5; i++)
    {
        verts.push_back(vec2(-.9, y));
        verts.push_back(vec2(.9, y));
        colours.push_back(vec3(0.f));
        colours.push_back(vec3(0.f));
        y -= 0.1;
    }


    // only one object selected
    if (selected1 != -1 && selected2 == -1)
    {
        Ray tempRay(ray->origin, ray->direction);
        objectVec[selected1]->getVolume(&tempRay);
        addToVerts(&tempRay, 0, objectVec[selected1]->colour);
    }
    // 2 objects selected
    else if (selected1 != -1 && selected2 != -1)
    {
        Ray tempRay(ray->origin, ray->direction);
        objectVec[selected1]->getVolume(&tempRay);
        addToVerts(&tempRay, 0, objectVec[selected1]->colour);

        objectVec[selected2]->getVolume(&tempRay);
        addToVerts(&tempRay, 1, objectVec[selected2]->colour);

        Object *obj = &Union(objectVec[selected1], objectVec[selected2]);
        obj->getVolume(&tempRay);
        addToVerts(&tempRay, 2, obj->colour);

        obj = new Intersection(objectVec[selected1], objectVec[selected2]);
        obj->getVolume(&tempRay);
        addToVerts(&tempRay, 3, obj->colour);

        obj = new Difference(objectVec[selected1], objectVec[selected2]);
        obj->getVolume(&tempRay);
        addToVerts(&tempRay, 4, obj->colour);

        objectVec[selected2]->differenceB = false;
    }

    // line scaling
    if (verts.size() > 10)
    {
        // get the minimum and maximum t values
        int i = 10;
        float tmin = FLT_MAX, tmax = FLT_MIN;
        do {
            tmin = min(tmin, verts[i].x);
            tmax = max(tmax, verts[i].x);
            i++;
        } while (i < verts.size());

        float distance = abs(tmax - tmin);

        for (int i = 10; i < verts.size(); i++)
            verts[i].x = (2.f * (.85f * (((verts[i].x - tmin) / distance))) + .15f) - 1.f;
    }

    glBindVertexArray(vertexArray);

    glBindBuffer(GL_ARRAY_BUFFER, vertexBuffer);
    glBufferData(GL_ARRAY_BUFFER, sizeof(verts[0]) * verts.size(), &verts[0], GL_STATIC_DRAW);
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 0, NULL);
    glEnableVertexAttribArray(0);

    glBindBuffer(GL_ARRAY_BUFFER, colourBuffer);
    glBufferData(GL_ARRAY_BUFFER, sizeof(colours[0]) * colours.size(), &colours[0], GL_STATIC_DRAW);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, NULL);
    glEnableVertexAttribArray(1);
 

    glBindVertexArray(0);
}