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
struct MyRay
{
	vec3 origin, direction;
	MyRay() {}
};

struct MyLight
{
	vec3 point, colour, ambient;
	MyLight() {}
};

struct MySphere
{
	vec3 center, colour, specular;
	float radius, phong, reflect;
	MySphere() {}
};

struct MyPlane
{
	vec3 normal, point, colour, specular;
	float phong, reflect;
	MyPlane() {}
};

struct MyTriangle
{
	float a, b, c, d, e, f, phong, reflect;
	vec3 p1, p2, p3, colour, specular;
	MyTriangle() {}
};

struct vertice
{
	float x, y, z;
	vertice() {}
};

struct polygon
{
	vertice vertA, vertB, vertC;
	float texA, TexB, TexC;
	vertice normA, normB, normC;

};


//------------------------------------
// GIVEN OPEN GL FUNCTIONS

void ErrorCallback(int error, const char* description);
void QueryGLVersion();
bool CheckGLErrors();

//------------------------------------
//LIGHTING FUNCTIONS

vec3 getColour(MyRay& ray, vector<MySphere>& sphereVec, vector<MyTriangle>& triangleVec, vector<MyPlane>& planeVec, vector<MyLight>& lightVec, int recursive);
vec3 shading(vec3 colourIn, vec3 intersection, vec3 origin, vector<MyLight>& lightVec, vec3 n, float phong, vec3 specular);
bool checkShadow(vec3 intersection, vector<MySphere>& sphereVec, vector<MyPlane>& planeVec, vector<MyTriangle>& triangleVec, vector<MyLight>& lightVec);

//------------------------------------
//SCALAR FUNCTIONS
	// spheres
float getNearestSphereScalar(MyRay ray, vector<MySphere>& sphereVec, vec3 *colourVec, vec3 *normal, float *phong, vec3 *specular, float *reflect);
float getNearestSphereScalar(MyRay ray, vector<MySphere>& sphereVec);
	//triangles
float getNearestTriangleScalar(MyRay ray, vector<MyTriangle>& triangleVec, vec3 *colourVec, vec3 *normal, float *phong, vec3 *specular, float *reflect);
float getNearestTriangleScalar(MyRay ray, vector<MyTriangle>& triangleVec);
	//planes
float getNearestPlaneScalar(MyRay ray, vector<MyPlane>& planeVec, vec3 *colourVec, vec3 *normal, float *phong, vec3 *specular, float *reflect);

//------------------------------------
// UTILITY FUNCTIONS
float magnitude(vec3 ray);
void myNormalize(vec3 *ray);
void myNormalize(MyRay& ray);

//-------------------------------------
// FILE FUNCTIONS
void readFromFile(const string fileDir, vector<MyLight>&  lightVec, vector<MySphere>& sphereVec, vector<MyPlane>& planeVec, vector<MyTriangle>& triangleVec);
void getPointFromLine(string line, vec3 *point);