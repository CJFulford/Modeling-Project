#include "Header.h"
#include "ImageBuffer.h"

void QueryGLVersion();
bool CheckGLErrors();

int scene = 1;
vec3 camOrigin = vec3(0.0, 0.0, 0.0);

void KeyCallback(
	GLFWwindow* window, 
	int key, 
	int scancode, 
	int action, 
	int mods)
{
    if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS)
        glfwSetWindowShouldClose(window, GL_TRUE);
	if (key == GLFW_KEY_1 && action == GLFW_PRESS)
		scene = 1;
	if (key == GLFW_KEY_2 && action == GLFW_PRESS)
		scene = 2;
	if (key == GLFW_KEY_3 && action == GLFW_PRESS)
		scene = 3;
	if (key == GLFW_KEY_LEFT && action == GLFW_PRESS)
		camOrigin += vec3(-0.1, 0.0, 0.0);
	if (key == GLFW_KEY_RIGHT && action == GLFW_PRESS)
		camOrigin += vec3(0.1, 0.0, 0.0);
	if (key == GLFW_KEY_UP && action == GLFW_PRESS)
		camOrigin += vec3(0.0, 0.1, 0.0);
	if (key == GLFW_KEY_DOWN && action == GLFW_PRESS)
		camOrigin += vec3(0.0, -0.1, 0.0);
	if (key == GLFW_KEY_O && action == GLFW_PRESS)
		camOrigin -= vec3(0.0, 0.0, 0.1);
	if (key == GLFW_KEY_P && action == GLFW_PRESS)
		camOrigin += vec3(0.0, 0.0, 0.1);

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
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 5);
	glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
	glfwWindowHint(GLFW_DOUBLEBUFFER, TRUE);
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

	int oldScene = scene - 1;

	vector<MyLight>		lightVec;
	vector<MySphere>	sphereVec;
	vector<MyPlane>		planeVec;
	vector<MyTriangle>	triangleVec;
	readFromFile("scene" + to_string(scene) + ".txt", lightVec, sphereVec, planeVec, triangleVec);

	while (!glfwWindowShouldClose(window))
	{
		//if (oldScene != scene) {
			// start recording time for render timer
			cout << "Rendering... ";
			time_t startTime = time(NULL);
			vec3 colourVec;

			// += 2.f since the window goes from -1 top 1, which is a distance of 2
			int xcoord;
			//#pragma omp parallel for schedule(dynamic)
			for (xcoord = -WINDOW_WIDTH; xcoord < WINDOW_WIDTH; xcoord += 2)
			{
				for (int ycoord = -WINDOW_HEIGHT; ycoord < WINDOW_HEIGHT; ycoord += 2)
				{
					int recursive = MAX_RECURSIVE_RAYS;
					MyRay ray;
					if (ratio >= 1.f) ray.direction = vec3((float)(xcoord / WINDOW_WIDTH) * ratio, (float)(ycoord / WINDOW_HEIGHT), FOCAL_LENGTH);
					else ray.direction = vec3(xcoord, ycoord / ratio, FOCAL_LENGTH);
					ray.origin = camOrigin;
					myNormalize(ray);
					colourVec = getColour(ray, sphereVec, triangleVec, planeVec, lightVec, recursive);
					//imageBuffer.SetPixel((xcoord * HALF_WIDTH) + HALF_WIDTH,	(ycoord * HALF_HEIGHT) + HALF_HEIGHT,	colourVec);
					imageBuffer.SetPixel(xcoord / 2 + HALF_WIDTH,	ycoord / 2 + HALF_HEIGHT,	colourVec);
				}
			}

			imageBuffer.Render();

			// scene is rendered to the back buffer, so swap to front for display
			glfwSwapBuffers(window);
			time_t endTime = time(NULL);
			cout << difftime(endTime, startTime) << " seconds"<< endl;
			oldScene = scene;
		//}
		// sleep until next event before drawing again
		glfwPollEvents();
	}

	// clean up allocated resources before exit
	glfwDestroyWindow(window);
	glfwTerminate();

	cout << "Goodbye!" << endl;
	return 0;
}
