#include "Header.h"
#include <iostream>

using namespace std;

int selected1 = -1, selected2 = -1;
double mouse_old_x, mouse_old_y;

float   rotate_x = 0.0,
        rotate_y = 0.0,
        zoom = DEF_ZOOM,
        aspectRatio = (float)WINDOW_WIDTH / (float)WINDOW_HEIGHT;

bool    select1   = true, 
        scale     = false, 

        movement  = false,
        movementX = false,
        movementY = false,
        movementZ = false,

        rotation  = false,
        rotationX = false,
        rotationY = false,
        rotationZ = false;
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
    if (action == GLFW_RELEASE)
    {
        switch (key)
        {
        case(GLFW_KEY_X):
            if (movement)      movementX = false;
            else if (rotation) rotationX = false;
            break;
        case(GLFW_KEY_Y):
            if (movement)      movementY = false;
            else if (rotation) rotationY = false;
            break;
        case(GLFW_KEY_Z):
            if (movement)      movementZ = false;
            else if (rotation) rotationZ = false;
            break;
        }
    }
	else if (action == GLFW_PRESS)
	{
		switch (key)
		{
		case(GLFW_KEY_ESCAPE):
			glfwSetWindowShouldClose(window, GL_TRUE);
			break;

        // object deletion
        case (GLFW_KEY_DELETE):
            if (selected1 != -1 && selected2 == -1)
            {
                delete objectVec[selected1];
                objectVec.erase(objectVec.begin() + selected1);
                selected1 = -1;
                select1 = true;
                scale = false;
                movement = false;
                rotation = false;
            }
            break;


        // scale
        case (GLFW_KEY_S):
            if (selected1 != -1 && selected2 == -1)
            {
                scale = !scale;
                movement = false;
                rotation = false;
            }
            break;
        // movement toggle
        case (GLFW_KEY_G):
            if (selected1 != -1 && selected2 == -1)
            {
                movement = !movement;
                scale = false;
                rotation = false;
            }
            break;
        // rotation toggle
        case (GLFW_KEY_R):
            if (selected1 != -1 && selected2 == -1)
            {
                rotation = !rotation;
                scale = false;
                movement = false;
            }
            break;


        // movement and rotation axis selections
        case(GLFW_KEY_X):
            if (movement)      movementX = true;
            else if (rotation) rotationX = true;
            break;
        case(GLFW_KEY_Y):
            if (movement)      movementY = true;
            else if (rotation) rotationY = true;
            break;
        case(GLFW_KEY_Z):
            if (movement)      movementZ = true;
            else if (rotation) rotationZ = true;
            break;


        
        // sphere
        case (GLFW_KEY_1):
            objectVec.push_back(new Sphere());
            break;
        // cube
        case (GLFW_KEY_2):
            objectVec.push_back(new Cube());
            break;
        case(GLFW_KEY_3):
            objectVec.push_back(new Torus());
            break;
        case(GLFW_KEY_4):
            objectVec.push_back(new Cylinder());
            break;
        
        
        case(GLFW_KEY_I): // intersection
            if (selected1 != -1 && selected2 != -1 && selected1 != selected2)
            {
                objectVec.push_back(new Intersection(objectVec[selected1], objectVec[selected2]));
                objectVec[selected1]->selected = false;
                objectVec[selected2]->selected = false;

                // need to check otherwise selected 2 will add the wrong object
                if (selected1 > selected2)
                {
                    objectVec.erase(objectVec.begin() + selected1);
                    objectVec.erase(objectVec.begin() + selected2);
                }
                else
                {
                    objectVec.erase(objectVec.begin() + selected1);
                    objectVec.erase(objectVec.begin() + selected2 - 1);
                }
                selected1 = -1;
                selected2 = -1;
                select1 = true;
            }
            break;
        case(GLFW_KEY_U): // union
            if (selected1 != -1 && selected2 != -1 && selected1 != selected2)
            {
                objectVec.push_back(new Union(objectVec[selected1], objectVec[selected2]));
                objectVec[selected1]->selected = false;
                objectVec[selected2]->selected = false;

                // need to check otherwise selected 2 will add the wrong object
                if (selected1 > selected2)
                {
                    objectVec.erase(objectVec.begin() + selected1);
                    objectVec.erase(objectVec.begin() + selected2);
                }
                else
                {
                    objectVec.erase(objectVec.begin() + selected1);
                    objectVec.erase(objectVec.begin() + selected2 - 1);
                }
                selected1 = -1;
                selected2 = -1;
                select1 = true;
            }
            break;
        case(GLFW_KEY_D): // difference
            if (selected1 != -1 && selected2 != -1 && selected1 != selected2)
            {
                objectVec.push_back(new Difference(objectVec[selected1], objectVec[selected2]));
                objectVec[selected1]->selected = false;
                objectVec[selected2]->selected = false;

                // need to check otherwise selected 2 will add the wrong object
                if (selected1 > selected2)
                {
                    objectVec.erase(objectVec.begin() + selected1);
                    objectVec.erase(objectVec.begin() + selected2);
                }
                else
                {
                    objectVec.erase(objectVec.begin() + selected1);
                    objectVec.erase(objectVec.begin() + selected2 - 1);
                }
                selected1 = -1;
                selected2 = -1;
                select1 = true;
            }
            break;
        case(GLFW_KEY_B): // Break join
            if (selected1 != -1 && selected2 == -1)
            {
                objectVec[selected1]->breakBoolean(&objectVec, selected1);
                selected1 = -1;
                select1 = true;
            }
            break;
		
        
        
        default:
			break;
		}
	}
}

void mouseMotion(GLFWwindow* window, double x, double y)
{
    #define SCREEN_CONTROL_SCALE (2.f / WINDOW_HEIGHT)

    // scale, rotation, movement applications
    if (scale)
    {
        if (y - mouse_old_y < 0)
            objectVec[selected1]->scale(true);
        else
            objectVec[selected1]->scale(false);
    }
    else if (movementX || rotationX)
    {
        if (movement)
            objectVec[selected1]->move(glm::vec3((float)(x - mouse_old_x) * SCREEN_CONTROL_SCALE, 0.f, 0.f));
        else if (rotation)
            objectVec[selected1]->rotate(glm::vec3((float)(x - mouse_old_x) * SCREEN_CONTROL_SCALE, 0.f, 0.f));
    }
    else if (movementY || rotationY)
    {
        if (movement)
            objectVec[selected1]->move(glm::vec3(0.f, (float)(y - mouse_old_y) * SCREEN_CONTROL_SCALE, 0.f));
        else if (rotation)
            objectVec[selected1]->rotate(glm::vec3(0.f, (float)(y - mouse_old_y) * SCREEN_CONTROL_SCALE, 0.f));
    }
    else if (movementZ || rotationZ)
    {
        if (movement)
            objectVec[selected1]->move(glm::vec3(0.f, 0.f, (float)(y - mouse_old_y) * SCREEN_CONTROL_SCALE));
        else if (rotation)
            objectVec[selected1]->rotate(glm::vec3(0.f, 0.f, (float)(y - mouse_old_y) * SCREEN_CONTROL_SCALE));
    }

    // screen rotation
    else if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_1))
    {
        rotate_x += (float)((y - mouse_old_y)) * radToDeg;
        rotate_y += (float)((x - mouse_old_x)) * radToDeg;
    }
    mouse_old_x = x;
    mouse_old_y = y;
}

void mouseButtonCallback(GLFWwindow* window, int button, int action, int mods)
{
    scale = false;
    movement = false;
    movementX = false;
    movementY = false;
    movementZ = false;
    rotation = false;
    rotationX = false;
    rotationY = false;
    rotationZ = false;
    // selection
    if (button == GLFW_MOUSE_BUTTON_RIGHT && action == GLFW_PRESS)
    {
        #define rayl	-1
        #define rayr	1
        #define rayt	1
        #define rayb	-1


        double xpos, ypos;
        glfwGetCursorPos(window, &xpos, &ypos);

        float u = rayl + ((rayr - rayl)*((float)xpos + .5f)) / WINDOW_WIDTH;
        float v = rayb + ((rayt - rayb)*(WINDOW_HEIGHT - (float)ypos + .5f)) / WINDOW_HEIGHT;
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
        glm::vec3 garbageNorm(0.f);
        Object *obj = ray.getClosestScalar(&scalar, &garbageNorm);

        auto temp = std::find(objectVec.begin(), objectVec.end(), obj);
        int index = temp - objectVec.begin();
        if (temp != objectVec.end() && 
            objectVec[index]->selectable && 
            selected1 != index && 
            selected2 != index)
        {
            if (select1)
            {
                if (selected1 != -1 && selected1 != index)
                    objectVec[selected1]->deselect();
                selected1 = index;
                objectVec[selected1]->select();
                select1 = false;
            }
            else
            {
                if (selected2 != -1 && selected2 != index)
                    objectVec[selected2]->deselect();
                selected2 = index;
                objectVec[selected2]->selected = true;
                select1 = true;
            }
        }
        else
        {
            if (selected1 != -1) objectVec[selected1]->deselect();
            if (selected2 != -1) objectVec[selected2]->deselect();
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

// %%%%%%%%%%%%%%% window creation
GLFWwindow* generateWindow()
{
	if (!glfwInit()) 
    {
		std::cout << "ERROR: GLFW failed to initilize, TERMINATING" << std::endl;
		return NULL;
	}

	// attempt to create a window with an OpenGL 4.4 core profile context
	GLFWwindow *window = glfwCreateWindow(WINDOW_WIDTH, WINDOW_HEIGHT, "Implicit Surfaces", 0, 0);
	if (!window) 
    {
		std::cout << "Program failed to create GLFW window, TERMINATING" << std::endl;
		glfwTerminate();
		return NULL;
	}

	// set Callbacks
    glfwSetErrorCallback(ErrorCallback);
    glfwSetKeyCallback(window, keyCallback);
    glfwSetCursorPosCallback(window, mouseMotion);
    glfwSetScrollCallback(window, scrollCallback);
    glfwSetMouseButtonCallback(window, mouseButtonCallback);
    glfwMakeContextCurrent(window);

	// load glad
	gladLoadGL();

	// query and print out information about our OpenGL environment
    printOpenGLVersion(GL_MAJOR_VERSION, GL_MINOR_VERSION, GL_SHADING_LANGUAGE_VERSION);
    
    return window;
}
