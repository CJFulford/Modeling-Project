#pragma once
#include <glm\glm.hpp>
#include <vector>

struct Object
{
	glm::vec3 colour;
	float phong;
};

struct Sphere : Object
{
	Sphere(glm::vec3 center, float radius, glm::vec3 col, float ph) : 
		center(center), radius(radius)
	{
		colour = col;
		phong = ph;
	}
	glm::vec3 center;
	float radius;
};

struct Triangle : Object
{
	Triangle(glm::vec3 p1, glm::vec3 p2, glm::vec3 p3, glm::vec3 col, float ph) : 
		p1(p1), p2(p2), p3(p3)
	{
		a = p1.x - p2.x;
		b = p1.y - p2.y;
		c = p1.z - p2.z;

		d = p1.x - p3.x;
		e = p1.y - p3.y;
		f = p1.z - p3.z;

		normal = normalize(cross(p3 - p2, p1 - p2));
		colour = col;
		phong = ph;
	}
	glm::vec3 p1, p2, p3, normal;
	float a, b, c, d, e, f;
};

struct Torus : Object
{
	Torus(glm::vec3 center, float mainRadius, float subRadius, glm::vec3 col, float ph) : 
		center(center), mainRadius(mainRadius), subRadius(subRadius) 
	{
		colour = col;
		phong = ph;
	}
	glm::vec3 center, colour;
	float mainRadius, subRadius, phong;
};

struct Volume
{
	Volume(glm::vec3 enterence, glm::vec3 exit) :
		enterence(enterence), exit(exit) {}
	glm::vec3 enterence, exit;
	Object object;
};

struct Ray
{
	Ray(glm::vec3 orig, glm::vec3 dir) :
		origin(orig), direction(dir) {}
	glm::vec3 origin, direction;
	std::vector<Volume> volumes;
};