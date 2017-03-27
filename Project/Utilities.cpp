#include "Header.h"
#include <iostream>

using namespace std;

double mouse_old_x, mouse_old_y;

float   rotate_x = 0.0,
        rotate_y = 0.0,
        zoom = DEF_ZOOM,
        aspectRatio = (float)WINDOW_WIDTH / (float)WINDOW_HEIGHT;

glm::vec3 rotation = DEF_ROTATION;

// Error Checking
void ErrorCallback(int error, const char* description)
{
	cout << "GLFW ERROR " << error << ":" << endl;
	cout << description << endl;
}

bool CheckGLErrors()
{
    bool error = false;
    for (GLenum flag = glGetError(); flag != GL_NO_ERROR; flag = glGetError())
    {
        cout << "OpenGL ERROR:  ";
        switch (flag) {
        case GL_INVALID_ENUM:
            cout << "GL_INVALID_ENUM" << endl; break;
        case GL_INVALID_VALUE:
            cout << "GL_INVALID_VALUE" << endl; break;
        case GL_INVALID_OPERATION:
            cout << "GL_INVALID_OPERATION" << endl; break;
        case GL_INVALID_FRAMEBUFFER_OPERATION:
            cout << "GL_INVALID_FRAMEBUFFER_OPERATION" << endl; break;
        case GL_OUT_OF_MEMORY:
            cout << "GL_OUT_OF_MEMORY" << endl; break;
        default:
            cout << "[unknown error code]" << endl;
        }
        error = true;
    }
    return error;
}


// For User Information
void printOpenGLVersion(GLenum majorVer, GLenum minorVer, GLenum langVer)
{
    GLint major, minor;
    glGetIntegerv(majorVer, &major);
    glGetIntegerv(minorVer, &minor);
    printf("OpenGL on %s %s\n", glGetString(GL_VENDOR), glGetString(GL_RENDERER));
    printf("OpenGL version supported %s\n", glGetString(GL_VERSION));
    printf("GLSL version supported %s\n", glGetString(langVer));
    printf("GL version major, minor: %i.%i\n", major, minor);
}


// controls
void keyCallback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
	if (action == GLFW_PRESS)
	{
		switch (key)
		{
		case(GLFW_KEY_ESCAPE):
			glfwSetWindowShouldClose(window, GL_TRUE);
			break;
		default:
			break;
		}
	}
}

void mouseMotion(GLFWwindow* window, double x, double y)
{
    if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_1))
    {
        rotate_x += (float)((y - mouse_old_y)) * radToDeg;
        rotate_y += (float)((x - mouse_old_x)) * radToDeg;
    }
    mouse_old_x = x;
    mouse_old_y = y;

}

void scrollCallback(GLFWwindow* window, double xoffset, double yoffset)
{
    if (yoffset < 0)
        zoom += ZOOM_SENS;
    else if (yoffset > 0)
        zoom = max(zoom - ZOOM_SENS, MAX_ZOOM);
}