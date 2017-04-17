#include "Shader.h"
#include <iostream>
#include <fstream>

unsigned long getFileLength(std::ifstream& file)
{
    if (!file.good()) return 0;

    file.seekg(0, std::ios::end);
    unsigned long len = (unsigned long)file.tellg();
    file.seekg(std::ios::beg);

    return len;
}

GLchar* loadshader(std::string filename)
{
    std::ifstream file;
    file.open(filename.c_str(), std::ios::in);
    if (!file)
    {
        std::cout << "FILE " << filename.c_str() << " NOT FOUND" << std::endl;
        return NULL;
    }
    unsigned long len = getFileLength(file);

    if (len == 0) return NULL;

    GLchar* ShaderSource = 0;

    ShaderSource = new char[len + 1];

    if (ShaderSource == 0) return NULL;

    ShaderSource[len] = 0;

    unsigned int i = 0;
    while (file.good())
    {
        ShaderSource[i] = file.get();
        if (!file.eof())
            i++;
    }

    ShaderSource[i] = 0;
    
    file.close();

    return ShaderSource;
}

void unloadshader(GLchar** ShaderSource)
{
    if (*ShaderSource != 0) delete[] * ShaderSource;
    *ShaderSource = 0;
}

void attachShader(GLuint &program, const char* fileName, GLuint shaderType)
{
    GLuint shader;
    const GLchar *shaderSource[] = { loadshader(fileName) };

    shader = glCreateShader(shaderType);
    glShaderSource(shader, 1, shaderSource, NULL);
    glCompileShader(shader);
    GLint status;

    glGetShaderiv(shader, GL_COMPILE_STATUS, &status);
    if (status == GL_FALSE)
    {
        GLint infoLogLength;
        glGetShaderiv(shader, GL_INFO_LOG_LENGTH, &infoLogLength);

        GLchar* strInfoLog = new GLchar[infoLogLength + 1];
        glGetShaderInfoLog(shader, infoLogLength, NULL, strInfoLog);

        std::cout << "\n\n" << fileName << std::endl;

        fprintf(stderr, "compilation error in shader vertex_shader: \n%s", strInfoLog);
        delete[] strInfoLog;
    }
    glAttachShader(program, shader);

    glDeleteShader(shader);
    unloadshader((GLchar**)shaderSource);
}

GLuint generateProgram(const char* vertexFilename, const char* fragmentFilename)
{
    GLuint program = glCreateProgram();
    attachShader(program, vertexFilename, GL_VERTEX_SHADER);
    attachShader(program, fragmentFilename, GL_FRAGMENT_SHADER);
    glLinkProgram(program);
    return program;
}

