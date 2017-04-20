#include "Header.h"
#include "texture.h"
#include <vector>
#include <string.h>

#pragma once
class Icon
{
public:
    Icon(std::string file);
	Icon(std::string file, glm::vec2 pos);

    std::vector<glm::vec2> positions;

	void render();
    void update(int positionIndex);

private:
    Texture texture;
    GLuint vertexBuffer;
    GLuint vertexArray;
    GLuint uvBuffer;
    GLuint textureID;
    GLuint program;

    const glm::vec2 defaultVerts[4] = {
        glm::vec2(-.03f, -.03f),
        glm::vec2(.03f, -.03f),
        glm::vec2(-.03f,  .03f),
        glm::vec2(.03f,  .03f)
    };
    glm::vec2 verts[4] = {
        glm::vec2(-.03f, -.03f),
        glm::vec2(.03f, -.03f),
        glm::vec2(-.03f,  .03f),
        glm::vec2(.03f,  .03f)
    };
    std::vector<glm::vec2> uvs;

    std::string filename;
    std::vector<unsigned char> image;
    unsigned int imageWidth;
    unsigned int imageHeight;

    void generateBuffer();
    void renderPosition();
    void loadImages();
};

