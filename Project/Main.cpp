#include "Header.h"
#include "ImageBuffer.h"
#include <glm/gtx/rotate_vector.hpp>

// definitions of variables used in ray tracing
#define rayl	-1
#define rayr	1
#define rayt	1
#define rayb	-1

glm::vec3 camOrigin = DEF_CAM_POS;

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
    glfwMakeContextCurrent(window);

	// load glad
	gladLoadGL();

	// query and print out information about our OpenGL environment
    printOpenGLVersion(GL_MAJOR_VERSION, GL_MINOR_VERSION, GL_SHADING_LANGUAGE_VERSION);

	ImageBuffer imageBuffer;
	imageBuffer.Initialize();

    std::vector<Sphere>		sphereVec;
    std::vector<Triangle>	triangleVec;
	readFromFile("scene1.txt", sphereVec, triangleVec);
	

	// variable initialization
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

				colourVec = getColour(ray, sphereVec, triangleVec);
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
	return 0;
}
