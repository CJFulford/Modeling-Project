#include "Header.h"

void ErrorCallback(
	int error, 
	const char* description)
{
	cout << "GLFW ERROR " << error << ":" << endl;
	cout << description << endl;
}

void QueryGLVersion()
{
	// query opengl version and renderer information
	string version = reinterpret_cast<const char *>(glGetString(GL_VERSION));
	string glslver = reinterpret_cast<const char *>(glGetString(GL_SHADING_LANGUAGE_VERSION));
	string renderer = reinterpret_cast<const char *>(glGetString(GL_RENDERER));

	cout << "OpenGL [ " << version << " ] "
		<< "with GLSL [ " << glslver << " ] "
		<< "on renderer [ " << renderer << " ]" << endl;
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