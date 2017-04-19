#include "CSGtree.h"

using namespace std;
using namespace glm;


CSGtree::CSGtree()
{
	vertexArray = 0;
	vertexBuffer = 0;
	program = generateProgram("CSGtree.vert", "CSGtree.frag");
	generateBuffer();
}


void CSGtree::generateBuffer()
{
	glGenVertexArrays(1, &vertexArray);
	glBindVertexArray(vertexArray);

	makeVerts();

	glGenBuffers(1, &vertexBuffer);
	glGenBuffers(1, &colourBuffer);

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
	{
		icons[i].position = verts[i];
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

void CSGtree::constructInfo(Object *obj, int level)		
{
	Object *temp;
	temp = dynamic_cast<Sphere*>(obj);
	if (temp != NULL) // then it is a sphere
	{
		info.push_back(vec2(level, 1));
		tempColours.push_back(temp->colour);
		Icon tempcon("icons/Sphere.png", vec2(0, 0));
		tempcon.loadImages();
		tempIcons.push_back(tempcon);
	}

	temp = dynamic_cast<Cube*>(obj);
	if (temp != NULL)
	{
		info.push_back(vec2(level, 2));
		tempColours.push_back(temp->colour);
		Icon tempcon("icons/Cube.png", vec2(0, 0));
		tempcon.loadImages();
		tempIcons.push_back(tempcon);
	}

	temp = dynamic_cast<Torus*>(obj);
	if (temp != NULL)
	{
		info.push_back(vec2(level, 3));
		tempColours.push_back(temp->colour);
		Icon tempcon("icons/Torus.png", vec2(0, 0));
		tempcon.loadImages();
		tempIcons.push_back(tempcon);
	}

	temp = dynamic_cast<Cylinder*>(obj);
	if (temp != NULL)
	{
		info.push_back(vec2(level, 4));
		tempColours.push_back(temp->colour);
		Icon tempcon("icons/Cylinder.png", vec2(0, 0));
		tempcon.loadImages();
		tempIcons.push_back(tempcon);
	}

	temp = dynamic_cast<Union*>(obj);
	if (temp != NULL)
	{
		info.push_back(vec2(level, 5));
		tempColours.push_back(temp->colour);
		constructInfo(dynamic_cast<Union*>(obj)->leftChild, level + 1);
		constructInfo(dynamic_cast<Union*>(obj)->rightChild, level + 1);
		Icon tempcon("icons/Union.png", vec2(0, 0));
		tempcon.loadImages();
		tempIcons.push_back(tempcon);
	}
		
	temp = dynamic_cast<Intersection*>(obj);
	if (temp != NULL)
	{
		info.push_back(vec2(level, 6));
		tempColours.push_back(temp->colour);
		constructInfo(dynamic_cast<Intersection*>(obj)->leftChild, level + 1);
		constructInfo(dynamic_cast<Intersection*>(obj)->rightChild, level + 1);
		Icon tempcon("icons/Intersection.png", vec2(0, 0));
		tempcon.loadImages();
		tempIcons.push_back(tempcon);
	}

	temp = dynamic_cast<Difference*>(obj);
	if (temp != NULL)
	{
		info.push_back(vec2(level, 7));
		tempColours.push_back(temp->colour);
		constructInfo(dynamic_cast<Difference*>(obj)->leftChild, level + 1);
		constructInfo(dynamic_cast<Difference*>(obj)->rightChild, level + 1);
		Icon tempcon("icons/Difference.png", vec2(0, 0));
		tempcon.loadImages();
		tempIcons.push_back(tempcon);
	}	
}

void CSGtree::makeVerts()
{
	// for nodes
	vector<int> levelTotalCount;

	for (vec2 vec : info)
	{
		if (vec.x > levelTotalCount.size())
			levelTotalCount.push_back(1);
		else
			levelTotalCount[vec.x-1]++;
	}

	vector<int> levelCount(levelTotalCount.size(),0);
	for (int i = 0; i < info.size(); i++)
	{
		levelCount[info[i].x - 1]++;
		verts.push_back(vec2(((float) levelCount[info[i].x-1]) / (float)(levelTotalCount[info[i].x-1] + 1), 1.f - (info[i].x / (levelTotalCount.size() + 1.f))));
		colours.push_back(tempColours[i]);
		icons.push_back(tempIcons[i]);
	}


	int leafctr = info.size() - 1;
	int parentctr = leafctr;

	while(leafctr > 0)
	{
		if (info[parentctr].x != info[leafctr].x - 1)
			parentctr--;
		else
		{
			verts.push_back(verts[parentctr]);
			verts.push_back(verts[leafctr]);
			colours.push_back(vec3(0.f));
			colours.push_back(vec3(0.f));
			leafctr--;
			parentctr = leafctr;
		}
	}


}