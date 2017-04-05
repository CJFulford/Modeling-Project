#include "Header.h"
#include <iostream>

using namespace std;

double mouse_old_x, mouse_old_y;

float   rotate_x = 0.0,
        rotate_y = 0.0,
        zoom = DEF_ZOOM,
        aspectRatio = (float)WINDOW_WIDTH / (float)WINDOW_HEIGHT;
bool select1 = true;
int selected1 = -1, selected2 = -1;

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

void mouseButtonCallback(GLFWwindow* window, int button, int action, int mods)
{
    if (button == GLFW_MOUSE_BUTTON_RIGHT && action == GLFW_PRESS)
    {
        #define rayl	-1
        #define rayr	1
        #define rayt	1
        #define rayb	-1


        double xpos, ypos;
        glfwGetCursorPos(window, &xpos, &ypos);

        float u = rayl + ((rayr - rayl)*(xpos + .5f)) / WINDOW_WIDTH;
        float v = rayb + ((rayt - rayb)*(WINDOW_HEIGHT - ypos + .5f)) / WINDOW_HEIGHT;
        float w = -(rayr / (float)tan(FOV / 2));

        // construct the ray
        // rotate the direction along the x axis then the y axis
        Ray ray(rotateY(    rotateX(camOrigin, rotate_x),     rotate_y) * zoom,
            glm::normalize(rotateY(    rotateX(glm::vec3(u, v, w), rotate_x),     rotate_y)));

        for (unsigned int i = 0; i < objectVec.size(); i++)
        {
            Object *object = (objectVec)[i];
            object->getVolume(&ray);
        }

        float scalar = 0;
        Object *obj = ray.getClosestScalar(&scalar);

        auto temp = std::find(objectVec.begin(), objectVec.end(), obj);
        if (temp != objectVec.end())
        {
            int index = temp - objectVec.begin();
            if (select1 && selected1 != index && selected2 != index)
            {
                if (selected1 != -1 && selected1 != index)
                    objectVec[selected1]->selected = false;
                selected1 = index;
                objectVec[selected1]->selected = true;
                std::cout << "\nsel1:" << selected1 << std::endl;
                select1 = false;
            }
            else if (selected1 != index && selected2 != index)
            {
                if (selected2 != -1 && selected2 != index)
                    objectVec[selected2]->selected = false;
                selected2 = index;
                objectVec[selected2]->selected = true;
                std::cout << "\nsel2:" << selected2 << std::endl;
                select1 = true;
            }
        }
        else
        {
            if (selected1 != -1) objectVec[selected1]->selected = false;
            if (selected2 != -1) objectVec[selected2]->selected = false;
            selected1 = -1;
            selected2 = -1;
            select1 = true;
        }
    }
}

void scrollCallback(GLFWwindow* window, double xoffset, double yoffset)
{
    if (yoffset < 0)
        zoom += ZOOM_SENS;
    else if (yoffset > 0)
        zoom = max(zoom - ZOOM_SENS, MAX_ZOOM);
}