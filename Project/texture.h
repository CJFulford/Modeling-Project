#include "Header.h"
#include <vector>
#include <string>

class Texture
{
public:
    Texture();

	 GLuint generateTexture(std::vector<unsigned char>& image, unsigned int width, unsigned int height);
	 void bindTexture(GLuint _program, GLuint _textureID, std::string varName);
	 void unbindTexture();	 
};
