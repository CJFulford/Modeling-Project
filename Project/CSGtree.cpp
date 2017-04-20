#include "CSGtree.h"

using namespace std;
using namespace glm;


CSGtree::CSGtree()
{
	vertexArray = 0;
	vertexBuffer = 0;
    colourBuffer = 0;
	program = generateProgram("CSGtree.vert", "CSGtree.frag");
	generateBuffer();
}


void CSGtree::generateBuffer()
{

	makeVerts();

    glBindVertexArray(vertexArray);
	glGenBuffers(1, &vertexBuffer);
	glGenBuffers(1, &colourBuffer);

    glGenVertexArrays(1, &vertexArray);

    glBindBuffer(GL_ARRAY_BUFFER, vertexBuffer);
    glBufferData(GL_ARRAY_BUFFER, sizeof(verts[0]) * verts.size(), &verts[0], GL_STATIC_DRAW);
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 0, NULL);
    glEnableVertexAttribArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, 0);

    glBindBuffer(GL_ARRAY_BUFFER, colourBuffer);
    glBufferData(GL_ARRAY_BUFFER, sizeof(colours[0]) * colours.size(), &colours[0], GL_STATIC_DRAW);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, NULL);
    glEnableVertexAttribArray(1);
    glBindBuffer(GL_ARRAY_BUFFER, 0);

    glBindVertexArray(0);

}

void CSGtree::render()
{
	glUseProgram(program);
	glBindVertexArray(vertexArray);

	//draw connecting lines
	glLineWidth(5);
	glDrawArrays(GL_LINES, info.size(), verts.size() - info.size());

	//draw points
	glPointSize(30.f);
	glDrawArrays(GL_POINTS, 0, info.size());

	//render the icon s
	for (int i = 0; i < icons.size(); i++)
	{
		icons[i].update();
		icons[i].render();
	}
	
	glBindVertexArray(0);
	glUseProgram(0);
}

void CSGtree::update()
{
	makeVerts();

	for (int i = 0; i < icons.size(); i++)
		icons[i].position = verts[i];

    // since the size of the vertex and colour buffer can increase, we need to recreate the buffer. 
    glDeleteBuffers(1, &vertexBuffer);
    glDeleteBuffers(1, &colourBuffer);

    glBindVertexArray(vertexArray);

    glGenBuffers(1, &vertexBuffer);
	glBindBuffer(GL_ARRAY_BUFFER, vertexBuffer);
	glBufferData(GL_ARRAY_BUFFER, sizeof(verts[0]) * verts.size(), &verts[0], GL_STATIC_DRAW);
	glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 0, NULL);
	glEnableVertexAttribArray(0);

    glGenBuffers(1, &colourBuffer);
	glBindBuffer(GL_ARRAY_BUFFER, colourBuffer);
	glBufferData(GL_ARRAY_BUFFER, sizeof(colours[0]) * colours.size(), &colours[0], GL_STATIC_DRAW);
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, NULL);
	glEnableVertexAttribArray(1);

	glBindVertexArray(0);
}

void CSGtree::constructInfo(Object *obj, int level)		
{
    switch (obj->objectID)
    {
    case(1):
    {
        info.push_back(vec2(level, 1));
        colours.push_back(obj->colour);
        Icon tempcon("icons/Sphere.png", vec2(0, 0));
        tempcon.loadImages();
        icons.push_back(tempcon);
        break;
    }
    case(2):
    {
        info.push_back(vec2(level, 2));
        colours.push_back(obj->colour);
        Icon tempcon("icons/Cube.png", vec2(0, 0));
        tempcon.loadImages();
        icons.push_back(tempcon);
        break;
    }
    case(3):
    {
        info.push_back(vec2(level, 3));
        colours.push_back(obj->colour);
        Icon tempcon("icons/Torus.png", vec2(0, 0));
        tempcon.loadImages();
        icons.push_back(tempcon);
        break;
    }
    case(4):
    {
        info.push_back(vec2(level, 4));
        colours.push_back(obj->colour);
        Icon tempcon("icons/Cylinder.png", vec2(0, 0));
        tempcon.loadImages();
        icons.push_back(tempcon);
        break;
    }
    case(5):
    {
        info.push_back(vec2(level, 5));
        colours.push_back(obj->colour);
        constructInfo(dynamic_cast<Union*>(obj)->leftChild, level + 1);
        constructInfo(dynamic_cast<Union*>(obj)->rightChild, level + 1);
        Icon tempcon("icons/Union.png", vec2(0, 0));
        tempcon.loadImages();
        icons.push_back(tempcon);
        break;
    }
    case(6):
    {
        info.push_back(vec2(level, 6));
        colours.push_back(obj->colour);
        constructInfo(dynamic_cast<Intersection*>(obj)->leftChild, level + 1);
        constructInfo(dynamic_cast<Intersection*>(obj)->rightChild, level + 1);
        Icon tempcon("icons/Intersection.png", vec2(0, 0));
        tempcon.loadImages();
        icons.push_back(tempcon);
        break;
    }
    case(7):
    {
        info.push_back(vec2(level, 7));
        colours.push_back(obj->colour);
        constructInfo(dynamic_cast<Difference*>(obj)->leftChild, level + 1);
        constructInfo(dynamic_cast<Difference*>(obj)->rightChild, level + 1);
        Icon tempcon("icons/Difference.png", vec2(0, 0));
        tempcon.loadImages();
        icons.push_back(tempcon);
        break;
    }
    default:
        break;
    }
}

void CSGtree::makeVerts()
{
	// for nodes
	vector<int> levelTotalCount;

    // count the number of objects in each level. index = level. value@index = # in level
	for (vec2 obj : info)
	{
        // indexing starts at 0, levels start at 1
		if (obj.x > levelTotalCount.size())
			levelTotalCount.push_back(1);
		else
			levelTotalCount[obj.x-1]++;
	}

    // generate the trees' nodes' positions
	vector<int> levelCount(levelTotalCount.size(),0);
	for (int i = 0; i < info.size(); i++)
	{
		levelCount[info[i].x - 1]++;    // indicate we have visited another object at CSG level x

        // based on the level and levelCound, generate the CSG x and y position
        float xcoord = (levelCount[info[i].x - 1]) / (levelTotalCount[info[i].x - 1] + 1.f);
        float ycoord = 1.f - (info[i].x / (levelTotalCount.size() + 1.f));

		verts.push_back(vec2(xcoord, ycoord));
	}


    // traverse the tree, generating the lines between the nodes
    // from each noe in the tree, traverse the tree in a reverse-depth first search and draw a line between that node and its parent
    // which will be the first node with a level of nodelevel-1. all lines are black
    for (int i = info.size() - 1; i >= 0; i--)
    {
        for (int j = i - 1; j < info.size(); j--)
        {
            if (info[j].x == info[i].x - 1)
            {
                verts.push_back(verts[j]);
                verts.push_back(verts[i]);
                colours.push_back(vec3(0.f));
                colours.push_back(vec3(0.f));
                break;
            }
        }
    }
}