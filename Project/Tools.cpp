#include "Tools.h"
#include <stdlib.h>

void printVec(glm::vec2 v)
{
	std::cout << v.x << "\t" << v.y << std::endl;
}
void printVec(glm::vec3 v)
{
	std::cout << v.x << "\t" << v.y << "\t" << v.z << std::endl;
}
void printVecArray(glm::vec2 *v, int size)
{
	for (int i = 0; i < size; i++)
		printVec(v[i]);
}
void printVecArray(glm::vec3 *v, int size)
{
	for (int i = 0; i < size; i++)
		printVec(v[i]);
}
void printVecVector(std::vector<glm::vec2> v)
{
	for (unsigned int i = 0; i < v.size(); i++)
		printVec(v[i]);
}
void printVecVector(std::vector<glm::vec3> v)
{
	for (unsigned int i = 0; i < v.size(); i++)
		printVec(v[i]);
}
void printFloatVector(std::vector<float> v)
{
    for (unsigned int i = 0; i < v.size(); i++)
        std::printf("%f\n", v[i]);
}

//Rodrigues' rotation formula
glm::vec3 rodriguesRotate(glm::vec3 vector, glm::vec3 axis, float angle)
{
	return glm::vec3((vector * (float)cos(angle)) +
		(glm::cross(axis, vector) * (float)sin(angle)) +
		(axis * glm::dot(axis, vector) * (1.f - (float)cos(angle))));
}

float randFloat(float startRange, float endRange)
{
     return startRange + (float)(rand()) / ((float)(RAND_MAX / (endRange - startRange)));
}

// generate random vector in range 0 - 1
glm::vec3 generateRandomVector()
{
    return glm::vec3(randFloat(0, 1), randFloat(0, 1), randFloat(0, 1));
}
// generate random vector
glm::vec3 generateRandomVector(float startRange, float endRange)
{
    return glm::vec3(randFloat(startRange, endRange), randFloat(startRange, endRange), randFloat(startRange, endRange));
}

// generate random vector
glm::vec3 generateRandomVector(float xStartRange, float xEndRange, float yStartRange, float yEndRange, float zStartRange, float zEndRange)
{
    return glm::vec3(randFloat(xStartRange, xEndRange), randFloat(yStartRange, yEndRange), randFloat(zStartRange, zEndRange));
}
