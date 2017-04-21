#include "Header.h"
#include <iostream>

using namespace std;

#define rayl	-1
#define rayr	1
#define rayt	1
#define rayb	-1

int selected1 = -1, selected2 = -1;
double mouse_old_x, mouse_old_y;
Ray tlistRay(DEF_CAM_POS, normalize(glm::vec3(0.f) - DEF_CAM_POS));

float   
    rotate_x = 0.0,
    rotate_y = 0.0,
    trotate_x = 0.0,
    trotate_y = 0.0,
    zoom = DEF_ZOOM,
    aspectRatio = (float)RENDER_WINDOW_WIDTH / (float)RENDER_WINDOW_HEIGHT;

bool
    select1 = true,
    scale = false,

    movement = false,
    rotation = false,

    alterX = false,
    alterY = false,
    alterZ = false;


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


// controls
void keyCallback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
    if (action == GLFW_RELEASE)
    {
        switch (key)
        {
        case(GLFW_KEY_X):
            if (movement | rotation) alterX = false;
            break;
        case(GLFW_KEY_Y):
            if (movement | rotation) alterY = false;
            break;
        case(GLFW_KEY_Z):
            if (movement | rotation) alterZ = false;
            break;
        case(GLFW_KEY_S):
            scale = false;
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
            if (selected1 != -1)
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
            if (movement | rotation) alterX = true;
            break;
        case(GLFW_KEY_Y):
            if (movement | rotation) alterY = true;
            break;
        case(GLFW_KEY_Z):
            if (movement | rotation) alterZ = true;
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
		
        case(GLFW_KEY_T):
        {
            double xpos, ypos;
            glfwGetCursorPos(window, &xpos, &ypos);

			xpos = 250.0;
			ypos = 250.0;

            float u = rayl + ((rayr - rayl)*((float)xpos + .5f)) / RENDER_WINDOW_WIDTH;
            float v = rayb + ((rayt - rayb)*(RENDER_WINDOW_HEIGHT - (float)ypos + .5f)) / RENDER_WINDOW_HEIGHT;
			//float u = HALF_RENDER_WIDTH;
			//float v = HALF_RENDER_HEIGHT;
            float w = -(rayr / (float)tan(FOV / 2));

            // construct the ray
            // rotate the direction along the x axis then the y axis
            tlistRay = Ray(rotateY(rotateX(camOrigin, rotate_x), rotate_y) * zoom,
                glm::normalize(rotateY(rotateX(glm::vec3(u, v, w), rotate_x), rotate_y)));

            trotate_x = rotate_x;
            trotate_y = rotate_y;

            break;
        }
        default:
			break;
		}
	}
}

void mouseMotion(GLFWwindow* window, double x, double y)
{
    #define SCREEN_CONTROL_SCALE (2.f / WHOLE_HEIGHT)

    // scale, rotation, movement applications
    if (scale)
    {
            objectVec[selected1]->scale(-(y - mouse_old_y) * SCREEN_CONTROL_SCALE);
            if (selected2 != -1)
                objectVec[selected2]->scale(-(y - mouse_old_y) * SCREEN_CONTROL_SCALE);
    }
    else if (alterX)
    {
		if (movement)
		{
			objectVec[selected1]->move(glm::vec3((float)(x - mouse_old_x) * SCREEN_CONTROL_SCALE, 0.f, 0.f));
			if(selected2 != -1)
				objectVec[selected2]->move(glm::vec3((float)(x - mouse_old_x) * SCREEN_CONTROL_SCALE, 0.f, 0.f));
		}  
        else if (rotation)
            objectVec[selected1]->rotate(glm::vec3((float)(x - mouse_old_x) * SCREEN_CONTROL_SCALE, 0.f, 0.f));
    }
    else if (alterY)
    {
		if (movement)
		{
			objectVec[selected1]->move(glm::vec3(0.f, (float)(y - mouse_old_y) * -SCREEN_CONTROL_SCALE, 0.f));
			if (selected2 != -1)
				objectVec[selected2]->move(glm::vec3(0.f, (float)(y - mouse_old_y) * -SCREEN_CONTROL_SCALE, 0.f));
		}
        else if (rotation)
            objectVec[selected1]->rotate(glm::vec3(0.f, (float)(y - mouse_old_y) * -SCREEN_CONTROL_SCALE, 0.f));
    }
    else if (alterZ)
    {
		if (movement)
		{
			objectVec[selected1]->move(glm::vec3(0.f, 0.f, (float)(y - mouse_old_y) * SCREEN_CONTROL_SCALE));
			if(selected2 != -1)
				objectVec[selected2]->move(glm::vec3(0.f, 0.f, (float)(y - mouse_old_y) * SCREEN_CONTROL_SCALE));
		}
        else if (rotation)
            objectVec[selected1]->rotate(glm::vec3(0.f, 0.f, (float)(y - mouse_old_y) * SCREEN_CONTROL_SCALE));
    }

    // screen rotation
    else if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_1))
    {
        if (x < WHOLE_WIDTH && y < WHOLE_HEIGHT)                  //moving around the Rendering
        {
            rotate_x += (float)((y - mouse_old_y)) * radToDeg;
            rotate_y += (float)((x - mouse_old_x)) * radToDeg;
        }
    }
    mouse_old_x = x;
    mouse_old_y = y;
}

void mouseButtonCallback(GLFWwindow* window, int button, int action, int mods)
{
    scale = false;
    movement = false;
    rotation = false;
    alterX = false;
    alterY = false;
    alterZ = false;
    // selection
    if (button == GLFW_MOUSE_BUTTON_RIGHT && action == GLFW_PRESS)
    {
        double xpos, ypos;
        glfwGetCursorPos(window, &xpos, &ypos);

        float u = rayl + ((rayr - rayl)*((float)xpos + .5f)) / RENDER_WINDOW_WIDTH;
        float v = rayb + ((rayt - rayb)*(RENDER_WINDOW_HEIGHT - (float)ypos + .5f)) / RENDER_WINDOW_HEIGHT;
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
void printControls()
{
    std::cout <<
        "1	create a sphere\n" <<
        "2	create a cube\n" <<
        "3	create a torus\n" <<
        "4 	create a cylinder\n" <<
        "U	perform the union operation on the 2 selected objects\n" <<
        "I	perform the intersection operation on the 2 selected objects\n" <<
        "D	perform the difference operation on the 2 selected objects*\n" <<
        "B	break the selected boolean object into its 2 children objects*\n" <<
        "DEL	delete the selected object\n" <<
        "G	enable object movement\n" <<
        "R	enable object rotation**\n" <<
        "S	enable object scaling**\n" <<
        "Hold X enable X-axis for movement/rotation***\n" <<
        "Hold Y enable Y-axis for movement/rotation***\n" <<
        "Hold Z enable Z-axis for movement/rotation***\n" <<
        "Mouse Controls\n" <<
        "\thold and drag the left mouse button to swivel the camera around the origin of the scene. movement, rotation, and scaling are disabled when this feature is used\n" <<
        "\tclick the right mouse button to select the object that the cursor is hovering over\n" <<
        "\twhen movement is enabled\n" <<
        "\t\tWhile holding X, move the mouse left and right to move the object along the X axis\n" <<
        "\t\tWhile holding Ys, move the mouse up and down and right to move the object along the Y axis\n" <<
        "\t\tWhile holding Z, move the mouse up and down to move the object along the Z axis\n" <<
        "\twhen Rotation is enabled\n" <<
        "\t\tWhile holding X, move the mouse left and right to rotate the object along the X axis\n" <<
        "\t\tWhile holding Ys, move the mouse up and down to rotate the object along the Y axis\n" <<
        "\t\tWhile holding Z, move the mouse up and down to rotate the object along the Z axis\n" <<
        "\twhen scaling is enabled move the mouse up and down respectively to increase and decrease the scaling of the object****\n" <<
        "\n" <<
        "*  	only a single object can be selected for this operation to work\n" <<
        "** 	only a single atomic object (Sphere, Cube, Torus, or Cylinder) can be selected for these functions to work\n" <<
        "***	holding X overrides holding Y which overrides holding Z\n" <<
        "****	scaling is currently uniform across an object. However, we have found that with enough */\n" << 
    std::endl;
}
GLFWwindow* generateWindow()
{
	if (!glfwInit()) 
    {
		std::cout << "ERROR: GLFW failed to initilize, TERMINATING" << std::endl;
		return NULL;
	}

    glfwWindowHint(GLFW_RESIZABLE, GLFW_FALSE);
    glfwWindowHint(GLFW_DOUBLEBUFFER, GLFW_TRUE);
	glfwWindowHint(GLFW_SAMPLES, 16);

	// attempt to create a window with an OpenGL 4.4 core profile context
	GLFWwindow *window = glfwCreateWindow(WHOLE_WIDTH, WHOLE_HEIGHT, "CSG Visualizer - Cody Fulford, Conlan Hanwell, and Brendan Petras", 0, 0);
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

	objectVec.push_back(new RayCylinder(new Ray(DEF_CAM_POS, normalize(ZERO_VECTOR - DEF_CAM_POS))));

    printControls();

    return window;
}

