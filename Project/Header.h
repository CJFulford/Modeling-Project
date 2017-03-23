#include <fstream>
#include <string>
#include <iterator>
#include <algorithm>
#include <vector>
#include <ctime>
#include <iostream>

#include <glad\glad.h>
#include <glcorearb.h>
#include <glm\glm.hpp>
#include <omp.h>

#define PI 3.1415
#define FOVdeg 60
#define FOV FOVdeg * PI / 180

#define WINDOW_WIDTH	500
#define WINDOW_HEIGHT	500
#define HALF_WIDTH		(WINDOW_WIDTH / 2)
#define HALF_HEIGHT		(WINDOW_HEIGHT / 2)
#define aspectRatio		((float)(WINDOW_WIDTH / WINDOW_HEIGHT));
#define FOCAL_LENGTH	-2.2f
#define RAY_RECURSIONS	0
#define FLOAT_ERROR		0.001f

#define BLACK			vec3(0.0, 0.0, 0.0)
#define WHITE			vec3(1.0, 1.0, 1.0);

using namespace std;
using namespace glm;

// shape structs
struct Ray
{
	vec3 origin, direction;
};

struct Light
{
	vec3 point, colour, ambient;
};

struct Sphere
{
	vec3 center, colour, specular;
	float radius, phong, reflect;
};

struct Triangle
{
	float a, b, c, d, e, f, phong, reflect;
	vec3 p1, p2, p3, colour, specular;
};


//------------------------------------
// GIVEN OPEN GL FUNCTIONS

void ErrorCallback(int error, const char* description);
void QueryGLVersion();
bool CheckGLErrors();

//------------------------------------
//LIGHTING FUNCTIONS

vec3 getColour(Ray& ray, vector<Sphere>& sphereVec, vector<Triangle>& triangleVec, vector<Light>& lightVec, int recursive);
vec3 shading(vec3 colourIn, vec3 intersection, vec3 origin, vector<Light>& lightVec, vec3 n, float phong, vec3 specular);

//------------------------------------
//SCALAR FUNCTIONS
	// spheres
float getNearestSphereScalar(Ray ray, vector<Sphere>& sphereVec, vec3 *colourVec, vec3 *normal, float *phong, vec3 *specular, float *reflect);
float getNearestSphereScalar(Ray ray, vector<Sphere>& sphereVec);
	//triangles
float getNearestTriangleScalar(Ray ray, vector<Triangle>& triangleVec, vec3 *colourVec, vec3 *normal, float *phong, vec3 *specular, float *reflect);
float getNearestTriangleScalar(Ray ray, vector<Triangle>& triangleVec);

//-------------------------------------
// FILE FUNCTIONS
void readFromFile(const string fileDir, vector<Light>&  lightVec, vector<Sphere>& sphereVec, vector<Triangle>& triangleVec);