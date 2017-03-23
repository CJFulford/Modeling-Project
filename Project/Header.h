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
using namespace std;
using namespace glm;

const int WINDOW_WIDTH = 500;
const int WINDOW_HEIGHT = WINDOW_WIDTH;
const int HALF_WIDTH = WINDOW_WIDTH / 2;
const int HALF_HEIGHT = WINDOW_HEIGHT / 2;
const float ratio = WINDOW_WIDTH / WINDOW_HEIGHT;
const float FOCAL_LENGTH = -2.2f;
const int MAX_RECURSIVE_RAYS = 0;
const float error = 0.001f;
const vec3 BLACK = vec3(0.0, 0.0, 0.0);
const vec3 WHITE = vec3(1.0, 1.0, 1.0);

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

struct Plane
{
	vec3 normal, point, colour, specular;
	float phong, reflect;
};

struct Triangle
{
	float a, b, c, d, e, f, phong, reflect;
	vec3 p1, p2, p3, colour, specular;
};

struct Vertex
{
	float x, y, z;
};


//------------------------------------
// GIVEN OPEN GL FUNCTIONS

void ErrorCallback(int error, const char* description);
void QueryGLVersion();
bool CheckGLErrors();

//------------------------------------
//LIGHTING FUNCTIONS

vec3 getColour(Ray& ray, vector<Sphere>& sphereVec, vector<Triangle>& triangleVec, vector<Plane>& planeVec, vector<Light>& lightVec, int recursive);
vec3 shading(vec3 colourIn, vec3 intersection, vec3 origin, vector<Light>& lightVec, vec3 n, float phong, vec3 specular);
bool checkShadow(vec3 intersection, vector<Sphere>& sphereVec, vector<Plane>& planeVec, vector<Triangle>& triangleVec, vector<Light>& lightVec);

//------------------------------------
//SCALAR FUNCTIONS
	// spheres
float getNearestSphereScalar(Ray ray, vector<Sphere>& sphereVec, vec3 *colourVec, vec3 *normal, float *phong, vec3 *specular, float *reflect);
float getNearestSphereScalar(Ray ray, vector<Sphere>& sphereVec);
	//triangles
float getNearestTriangleScalar(Ray ray, vector<Triangle>& triangleVec, vec3 *colourVec, vec3 *normal, float *phong, vec3 *specular, float *reflect);
float getNearestTriangleScalar(Ray ray, vector<Triangle>& triangleVec);
	//planes
float getNearestPlaneScalar(Ray ray, vector<Plane>& planeVec, vec3 *colourVec, vec3 *normal, float *phong, vec3 *specular, float *reflect);

//-------------------------------------
// FILE FUNCTIONS
void readFromFile(const string fileDir, vector<Light>&  lightVec, vector<Sphere>& sphereVec, vector<Plane>& planeVec, vector<Triangle>& triangleVec);
void getPointFromLine(string line, vec3 *point);