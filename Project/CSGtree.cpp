#include "CSGtree.h"
#include "Shader.h"

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

	//getLines(NULL);

	glGenBuffers(1, &vertexBuffer);
	glBindBuffer(GL_ARRAY_BUFFER, vertexBuffer);
	glBufferData(GL_ARRAY_BUFFER, sizeof(verts[0]) * verts.size(), &verts[0], GL_STATIC_DRAW);
	glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 0, NULL);
	glEnableVertexAttribArray(0);
	glBindVertexArray(0);
}

void CSGtree::render()
{
	glUseProgram(program);
	glBindVertexArray(vertexArray);

	//draw points
	glPointSize(2.f);
	glDrawArrays(GL_POINTS, 0, verts.size());

	//draw connecting lines
	//glLineWidth(3);
	//glDrawArrays(GL_LINES, 0, 10);                 

	glBindVertexArray(0);
	glUseProgram(0);

}

void CSGtree::update()
{
	//getLines(NULL);

	glBindVertexArray(vertexArray);

	glBindBuffer(GL_ARRAY_BUFFER, vertexBuffer);
	glBufferData(GL_ARRAY_BUFFER, sizeof(verts[0]) * verts.size(), &verts[0], GL_STATIC_DRAW);
	glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 0, NULL);
	glEnableVertexAttribArray(0);

	glBindVertexArray(0);
}

void CSGtree::constructInfo(Object *obj)
{
	obj
}