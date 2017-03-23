#include "Header.h"
#include "ImageBuffer.h"


#define rayl	-1
#define rayr	1
#define rayt	1
#define rayb	-1

vec3 camOrigin = vec3(0.0, 0.0, 0.0);

int main(
	int argc, 
	char *argv[])
{
	// initialize the GLFW windowing system
	if (!glfwInit()) {
		std::cout << "ERROR: GLFW failed to initilize, TERMINATING" << endl;
		return -1;
	}
	glfwSetErrorCallback(ErrorCallback);

	// attempt to create a window with an OpenGL 4.1 core profile context
	GLFWwindow *window;
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 4);
	glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
	window = glfwCreateWindow(WINDOW_WIDTH, WINDOW_HEIGHT, "Ray Tracer Window", 0, 0);
	if (!window) {
		std::cout << "Program failed to create GLFW window, TERMINATING" << endl;
		glfwTerminate();
		return -1;
	}

	// set keyboard callback function and make our context current (active)
	glfwSetKeyCallback(window, KeyCallback);
	glfwMakeContextCurrent(window);

	// load glad
	gladLoadGL();

	// query and print out information about our OpenGL environment
	QueryGLVersion();

	ImageBuffer imageBuffer;
	imageBuffer.Initialize();

	vector<Sphere>		sphereVec;
	vector<Triangle>	triangleVec;
	readFromFile("scene1.txt", sphereVec, triangleVec);
	

	// variable initialization
	double diff = 0;
	time_t startTime = 0, endTime = 0;
	vec3 colourVec = BLACK;
	int frames = 0;
	float	w = -(rayr / (float)tan(FOV / 2)),
			u = 0,
			v = 0;


	while (!glfwWindowShouldClose(window))
	{
		startTime = time(NULL);
		colourVec = BLACK;


		
		#pragma omp parallel for schedule(dynamic)
		for (int i = 0; i < WINDOW_WIDTH; i++)
		{
			for (int j = 0; j < WINDOW_HEIGHT; j++)
			{
				Ray ray;

				u = rayl + ((rayr - rayl)*(i + .5f)) / WINDOW_WIDTH;
				v = rayb + ((rayt - rayb)*(j + .5f)) / WINDOW_HEIGHT;

				ray.origin = camOrigin;
				ray.direction = normalize(vec3(u, v, w) - ray.origin);

				colourVec = getColour(ray, sphereVec, triangleVec);
				imageBuffer.SetPixel(i, j, colourVec);
			}
		}



		imageBuffer.Render();



		// scene is rendered to the back buffer, so swap to front for display
		glfwSwapBuffers(window);
		endTime = time(NULL);
		diff = difftime(endTime, startTime);
		frames++;
		if (diff >= 1)
		{
			std::cout << frames << std::endl;
			frames = 0;
		}




		glfwPollEvents();
	}

	// clean up allocated resources before exit
	glfwDestroyWindow(window);
	glfwTerminate();

	std::cout << "Goodbye!" << endl;
	return 0;
}
