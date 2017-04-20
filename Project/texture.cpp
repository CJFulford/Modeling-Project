#include "texture.h"
#include <cmath>
#include <iostream>

Texture::Texture() {}

GLuint Texture::generateTexture(std::vector<unsigned char>& image, unsigned int width, unsigned int height)
{
	GLuint textureID;

	glGenTextures(1, &textureID);
	glBindTexture(GL_TEXTURE_2D, textureID);
	glTexStorage2D(GL_TEXTURE_2D, 1, GL_RGBA8, width, height);
	glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, width, height, GL_RGBA, GL_UNSIGNED_BYTE, image.data());
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);

	return textureID;
}

void Texture::bindTexture(GLuint _program, GLuint _textureID, std::string varName)
{
	glActiveTexture(GL_TEXTURE0 + _textureID);
	glBindTexture(GL_TEXTURE_2D, _textureID);
	glUniform1i(glGetUniformLocation(_program, varName.c_str()), _textureID);
}

void Texture::unbindTexture() { glBindTexture(GL_TEXTURE_2D, 0); }
