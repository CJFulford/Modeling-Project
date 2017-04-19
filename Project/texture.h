#ifndef TEXTURE_H
#define TEXTURE_H
#include "Header.h"
#include <vector>
#include <string>

class Texture
{
public:
    Texture();

	 GLuint create1DTexture(std::vector<GLubyte>& rgbaValues);
	 GLuint create2DTexture(std::vector<unsigned char>& image, unsigned int width, unsigned int height);
	 
	 void bind1DTexture(GLuint _program, GLuint _textureID, std::string varName);
	 void bind2DTexture(GLuint _program, GLuint _textureID, std::string varName);
	 
	 void unbind1DTexture();
	 void unbind2DTexture();	 

};

#endif // TEXTURE_H
