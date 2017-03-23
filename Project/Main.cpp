#include "Header.h"
#include "ImageBuffer.h"

void QueryGLVersion();
bool CheckGLErrors();

vec3 camOrigin = vec3(0.0, 0.0, 0.0);

void KeyCallback(
	GLFWwindow* window, 
	int key, 
	int scancode, 
	int action, 
	int mods)
{
	if (action == GLFW_PRESS)
	{
		switch (key)
		{
		case(GLFW_KEY_ESCAPE):
			glfwSetWindowShouldClose(window, GL_TRUE);

		case (GLFW_KEY_LEFT):
			camOrigin -= vec3(-0.1, 0.0, 0.0);
			break;
		case (GLFW_KEY_RIGHT):
			camOrigin -= vec3(0.1, 0.0, 0.0);
			break;
		case (GLFW_KEY_UP):
			camOrigin -= vec3(0.0, 0.1, 0.0);
			break;
		case(GLFW_KEY_DOWN):
			camOrigin -= vec3(0.0, -0.1, 0.0);
			break;
		case(GLFW_KEY_O):
			camOrigin -= vec3(0.0, 0.0, 0.1);
			break;
		case(GLFW_KEY_P):
			camOrigin += vec3(0.0, 0.0, 0.1);
			break;

		default:
			break;
		}
	}
}

int main(
	int argc, 
	char *argv[])
{
	// initialize the GLFW windowing system
	if (!glfwInit()) {
		cout << "ERROR: GLFW failed to initilize, TERMINATING" << endl;
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
		cout << "Program failed to create GLFW window, TERMINATING" << endl;
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

	vector<Light>		lightVec;
	vector<Sphere>	sphereVec;
	vector<Plane>		planeVec;
	vector<Triangle>	triangleVec;
	readFromFile("scene1.txt", lightVec, sphereVec, planeVec, triangleVec);

	while (!glfwWindowShouldClose(window))
	{
		// start recording time for render timer
		cout << "Rendering... ";
		time_t startTime = time(NULL);
		vec3 colourVec;



		int l = -1, r = 1, t = 1, b = -1;
		float w = -(r / tan(FOV / 2));
		/*float	u = l + ((r - l)*(i + .5f)) / WINDOW_WIDTH,
		v = b + ((t - b)*(j + .5f)) / WINDOW_HEIGHT;*/

		// += 2.f since the window goes from -1 top 1, which is a distance of 2
		#pragma omp parallel for schedule(dynamic)
		for (int i = 0; i < WINDOW_WIDTH; i++)
		{
			for (int j = 0; j < WINDOW_HEIGHT; j++)
			{
				int recursive = MAX_RECURSIVE_RAYS;
				Ray ray;

				float	u = -1 + (2*i + 1.f) / WINDOW_WIDTH,
						v = -1 + (2*j + 1.f) / WINDOW_HEIGHT;

				ray.origin = camOrigin;
				ray.direction = normalize(vec3(u, v, w) - ray.origin);

				colourVec = getColour(ray, sphereVec, triangleVec, planeVec, lightVec, recursive);
				imageBuffer.SetPixel(i, j, colourVec);
			}
		}

		imageBuffer.Render();

		// scene is rendered to the back buffer, so swap to front for display
		glfwSwapBuffers(window);
		time_t endTime = time(NULL);
		cout << difftime(endTime, startTime) << " seconds"<< endl;
		// sleep until next event before drawing again
		glfwPollEvents();
	}

	// clean up allocated resources before exit
	glfwDestroyWindow(window);
	glfwTerminate();

	cout << "Goodbye!" << endl;
	return 0;
}
