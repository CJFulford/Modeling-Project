#include <fstream>
#include <string>
#include <iterator>
#include <algorithm>
#include <vector>
#include <ctime>
#include <iostream>

#include <glad\glad.h>
#include <GLFW\glfw3.h>
#include <glm\glm.hpp>
#include <omp.h>

// Mathematical values
#define PI				3.14159265359f
#define radToDeg        (PI / 180)
#define identity		glm::mat4(1.f)
#define FLOAT_ERROR		0.001f

// window info
#define WINDOW_WIDTH	500
#define WINDOW_HEIGHT	500
#define HALF_WIDTH		(WINDOW_WIDTH / 2)
#define HALF_HEIGHT		(WINDOW_HEIGHT / 2)

// basic colours
#define BLACK			glm::vec3(0.f, 0.f, 0.f)
#define WHITE			glm::vec3(1.f, 1.f, 1.f)
#define RED             glm::vec3(1.f, 0.f, 0.f)
#define GREEN           glm::vec3(0.f, 1.f, 0.f)
#define BLUE            glm::vec3(0.f, 0.f, 1.f)

// camera info
#define FOV_DEGREE  	45
#define FOV				FOV_DEGREE * PI / 180
#define FOCAL_LENGTH	-2.2f
#define DEF_CAM_POS     glm::vec3(0.f, 0.f, 4.f)
#define ZOOM_SENS       .1f
#define DEF_ZOOM        1.f
#define MAX_ZOOM        .001f
#define DEF_ROTATION    glm::vec3(0.f, 0.f, 0.f)

extern glm::vec3 camOrigin;
extern glm::vec3 rotation;
extern float zoom;
extern float rotate_x;
extern float rotate_y;

// shape structs
struct Ray
{
    Ray(glm::vec3 orig, glm::vec3 dir) : origin(orig), direction(dir) {}
	glm::vec3 origin, direction;
};

struct Sphere
{
    glm::vec3 center, colour, specular;
	float radius, phong, reflect;
};

struct Triangle
{
    glm::vec3 p1, p2, p3, colour, specular;
	float a, b, c, d, e, f, phong, reflect;
};

// ============================== Utilities.cpp
// Error Checking
void ErrorCallback(int error, const char* description);
bool CheckGLErrors();

// user info
void printOpenGLVersion(GLenum majorVer, GLenum minorVer, GLenum langVer);

// control Functions
void keyCallback(GLFWwindow* window, int key, int scancode, int action, int mods);
void scrollCallback (GLFWwindow* window, double xoffset, double yoffset);
void mouseMotion(GLFWwindow* window, double x, double y);
void printOpenGLVersion(GLenum majorVer, GLenum minorVer, GLenum langVer);

// ============================== Lighting.cpp
glm::vec3 getColour(Ray& ray, std::vector<Sphere>& sphereVec, std::vector<Triangle>& triangleVec);
glm::vec3 shading(glm::vec3 colourIn, glm::vec3 intersection, glm::vec3 origin, glm::vec3 n, float phong, glm::vec3 specular);

// ============================== Scalars.cpp
float getNearestSphereScalar(Ray ray, std::vector<Sphere>& sphereVec, glm::vec3 *colourVec, glm::vec3 *normal, float *phong, glm::vec3 *specular, float *reflect);
float getNearestTriangleScalar(Ray ray, std::vector<Triangle>& triangleVec, glm::vec3 *colourVec, glm::vec3 *normal, float *phong, glm::vec3 *specular, float *reflect);

// ============================== File.cpp
void readFromFile(const std::string fileDir, std::vector<Sphere>& sphereVec, std::vector<Triangle>& triangleVec);