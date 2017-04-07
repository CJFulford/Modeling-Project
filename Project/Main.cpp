#include "Header.h"
#include "ImageBuffer.h"

#include <iostream>
#include <ctime>
#include <omp.h>
#include <glm/gtx/rotate_vector.hpp>

// definitions of variables used in ray tracing


glm::vec3 camOrigin = DEF_CAM_POS;
std::vector<Object*> objectVec;

// default scene consists of 2 spheres and a pyramid
void addFloor(std::vector<Object*> *objectVec)
{
    objectVec->push_back(new Torus(
        glm::vec3(0.f, 0.f, 0.f),           // center
        .5f,                                // major radius
        .2f,                                // minor radius
        glm::vec3(0.f, 1.f, 0.f),           // colour
        50.f));                             // phong
    objectVec->push_back(new Triangle(
        glm::vec3(-2.f, -2.f, -2.f),        //p1
        glm::vec3(-2.f, -2.f, 2.f),         //p2
        glm::vec3( 2.f, -2.f, -2.f),        //p3
        glm::vec3(0.1f, 0.8f, 0.9f),        //colour
        50.f));                             //phong
    objectVec->push_back(new Triangle(
        glm::vec3( 2.f, -2.f, -2.f),
        glm::vec3(-2.f, -2.f,  2.f),
        glm::vec3( 2.f, -2.f,  2.f),
        glm::vec3(0.1f, 0.8f, 0.9f),
        50.f));
}

int main(int argc, char *argv[])
{
	// initialize the GLFW windowing system
	if (!glfwInit()) 
    {
		std::cout << "ERROR: GLFW failed to initilize, TERMINATING" << std::endl;
		return -1;
	}

	// attempt to create a window with an OpenGL 4.4 core profile context
	GLFWwindow *window = glfwCreateWindow(WINDOW_WIDTH, WINDOW_HEIGHT, "Implicit Surfaces", 0, 0);
	if (!window) 
    {
		std::cout << "Program failed to create GLFW window, TERMINATING" << std::endl;
		glfwTerminate();
		return -1;
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

	ImageBuffer imageBuffer;
	imageBuffer.Initialize();

    addFloor(&objectVec);
	

	// variable initialization
    #define rayl	-1
    #define rayr	1
    #define rayt	1
    #define rayb	-1
	time_t startTime = 0, endTime = 0;
    glm::vec3 colourVec = BLACK;
	int frames = 0;
	float	w = -(rayr / (float)tan(FOV / 2)),
			u = 0,
			v = 0;


	while (!glfwWindowShouldClose(window))
	{
        
        startTime = time(NULL);

        // rotateX, rotateY, zoom. Do here for single calc
        glm::vec3 origin =  rotateY(
                                rotateX(camOrigin, rotate_x)
                                , rotate_y) 
                            * zoom;

		#pragma omp parallel for schedule(dynamic)
		for (int i = 0; i < WINDOW_WIDTH; i++)
		{
			for (int j = 0; j < WINDOW_HEIGHT; j++)
			{

				u = rayl + ((rayr - rayl)*(i + .5f)) / WINDOW_WIDTH;
				v = rayb + ((rayt - rayb)*(j + .5f)) / WINDOW_HEIGHT;

                // construct the ray
                // rotate the direction along the x axis then the y axis
				Ray ray(origin, 
                        glm::normalize(
                            rotateY(
                                rotateX(glm::vec3(u, v, w), rotate_x)
                                , rotate_y)));

				colourVec = getColour(&ray, &objectVec);
				imageBuffer.SetPixel(i, j, colourVec);
			}
		}



		imageBuffer.Render();
        glfwSwapBuffers(window);



        // update and print frames per second
		endTime = time(NULL);
		frames++;
		if (difftime(endTime, startTime) >= 1)
		{
            printf("\rFPS: %i", frames);
			frames = 0;
		}




		glfwPollEvents();
	}

	// clean up allocated resources before exit
	glfwDestroyWindow(window);
	glfwTerminate();
    for (Object *object : objectVec)
    {
        delete object;
    }
	return 0;
}
