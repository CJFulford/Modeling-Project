#include "Objects.h"
#include <glad/glad.h>
#include <GLFW/glfw3.h>

// window info
#define RENDER_WINDOW_WIDTH	    500
#define RENDER_WINDOW_HEIGHT	500
#define HALF_RENDER_WIDTH		(RENDER_WINDOW_WIDTH / 2)
#define HALF_RENDER_HEIGHT		(RENDER_WINDOW_HEIGHT / 2)
#define WHOLE_WIDTH             1000
#define WHOLE_HEIGHT            700
#define HALF_WHOLE_WIDTH		(WHOLE_WIDTH / 2)
#define HALF_WHOLE_HEIGHT		(WHOLE_HEIGHT / 2)

// basic colours
#define BLACK			glm::vec3(0.f)
#define WHITE			glm::vec3(1.f)
#define RED             glm::vec3(1.f, 0.f, 0.f)
#define GREEN           glm::vec3(0.f, 1.f, 0.f)
#define BLUE            glm::vec3(0.f, 0.f, 1.f)
#define AMBIENT         glm::vec3(0.3f)

// camera info
#define FOV_DEGREE  	45
#define FOV				FOV_DEGREE * PI / 180
#define FOCAL_LENGTH	-2.2f
#define DEF_CAM_POS     glm::vec3(0.f, 0.f, 4.f)
#define ZOOM_SENS       .1f
#define DEF_ZOOM        1.5f
#define MAX_ZOOM        .001f
#define DEF_ROTATION    glm::vec3(0.f)
#define DEF_MOVEMENT	glm::vec3(0.f)

#define LIGHT_POS       glm::vec3(2.f, 2.f, 2.f)
#define BACKGROUND_COLOUR   glm::vec3(.5f);

extern std::vector<Object*> objectVec;
extern glm::vec3 camOrigin;
extern float zoom;
extern Ray tlistRay;
extern int selected1;
extern int selected2;


 
// ============================== Utilities.cpp
// Error Checking
GLFWwindow* generateWindow();

// ============================== Lighting.cpp
glm::vec3 getColour(Ray *ray, std::vector<Object*> *objectVec);
