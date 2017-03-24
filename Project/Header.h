#include <vector>

#include <glad\glad.h>
#include <GLFW\glfw3.h>
#include <glm\glm.hpp>

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
#define BLACK			glm::vec3(0.f)
#define WHITE			glm::vec3(1.f)
#define RED             glm::vec3(1.f, 0.f, 0.f)
#define GREEN           glm::vec3(0.f, 1.f, 0.f)
#define BLUE            glm::vec3(0.f, 0.f, 1.f)
#define AMBIENT         glm::vec3(0.2f)

// camera info
#define FOV_DEGREE  	45
#define FOV				FOV_DEGREE * PI / 180
#define FOCAL_LENGTH	-2.2f
#define DEF_CAM_POS     glm::vec3(0.f, 0.f, 4.f)
#define ZOOM_SENS       .1f
#define DEF_ZOOM        1.f
#define MAX_ZOOM        .001f
#define DEF_ROTATION    glm::vec3(0.f)

#define LIGHT_POS       glm::vec3(1.f, 2.f, 1.f)

extern glm::vec3 camOrigin;
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
    Sphere(glm::vec3 center, float radius, glm::vec3 colour, float phong) : center(center), colour(colour), radius(radius), phong(phong) {}
    glm::vec3 center, colour;
	float radius, phong;
};

struct Triangle
{
    Triangle(glm::vec3 p1, glm::vec3 p2, glm::vec3 p3, glm::vec3 colour, float phong) : p1(p1), p2(p2), p3(p3), colour(colour), phong(phong) 
    {
        a = p1.x - p2.x;
        b = p1.y - p2.y;
        c = p1.z - p2.z;

        d = p1.x - p3.x;
        e = p1.y - p3.y;
        f = p1.z - p3.z;

        normal = normalize(cross(p3 - p2, p1 - p2));
    }
    glm::vec3 p1, p2, p3, colour, normal;
	float phong, a, b, c, d, e, f;
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
glm::vec3 getColour(Ray *ray, std::vector<Sphere> *sphereVec, std::vector<Triangle> *triangleVec);
glm::vec3 Blinn_Phong(Ray *ray, float scalar, glm::vec3 colourIn, glm::vec3 norm, float phong);

// ============================== Scalars.cpp
float getSphereScalar(Ray *ray, std::vector<Sphere> *sphereVec, glm::vec3 *normal, unsigned int *index);
float getTriScalar(Ray *ray, std::vector<Triangle> *triangleVec, unsigned int *index);
