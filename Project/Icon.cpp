#include "Icon.h"
#include "lodepng.h"
#include "ShaderBuilder.h"

using namespace std;
using namespace glm;



Icon::Icon(std::string file)
{
    vertexArray = 0;
    vertexBuffer = 0;
    uvBuffer = 0;
    textureID = 0;
    filename = file;
    positions;
    // all of our images are 32x32
    imageHeight = 32;
    imageWidth = 32;
    // shader program
    program = generateProgram("Icon.vert", "Icon.frag");
    // uv coordinates
    uvs.push_back(vec2(0.f, 1.f));
    uvs.push_back(vec2(1.f, 1.f));
    uvs.push_back(vec2(0.f, 0.f));
    uvs.push_back(vec2(1.f, 0.f));
    loadImages();
}

Icon::Icon(string file, vec2 pos)
{
	vertexArray = 0;
	vertexBuffer = 0;
	uvBuffer = 0;
	textureID = 0;
	filename = file;
	positions.push_back(pos);
    // all of our images are 32x32
	imageHeight = 32;
	imageWidth = 32;
    // shader program
	program = generateProgram("Icon.vert", "Icon.frag");
    // uv coordinates
	uvs.push_back(vec2(0.f, 1.f));
	uvs.push_back(vec2(1.f, 1.f));
	uvs.push_back(vec2(0.f, 0.f));
	uvs.push_back(vec2(1.f, 0.f));
    loadImages();
}

void Icon::generateBuffer()
{
	glGenVertexArrays(1, &vertexArray);
	glBindVertexArray(vertexArray);

	glGenBuffers(1, &vertexBuffer);
    glBindBuffer(GL_ARRAY_BUFFER, vertexBuffer);
    glBufferData(GL_ARRAY_BUFFER, sizeof(verts[0]) * sizeof(verts), verts, GL_STATIC_DRAW);
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 0, NULL);
    glEnableVertexAttribArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, 0);

    glGenBuffers(1, &uvBuffer);
    glBindBuffer(GL_ARRAY_BUFFER, uvBuffer);
    glBufferData(GL_ARRAY_BUFFER, sizeof(vec2) * uvs.size(), &uvs[0], GL_STATIC_DRAW);
    glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 0, NULL);
    glEnableVertexAttribArray(1);
    glBindBuffer(GL_ARRAY_BUFFER, 0);

    glBindVertexArray(0);
}

void Icon::renderPosition()
{
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    glUseProgram(program);
    glBindVertexArray(vertexArray);
    texture.bindTexture(program, textureID, std::string("tex"));

    //draw 2 triangles 
    glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);

    glBindVertexArray(0);
    texture.unbindTexture();
    glUseProgram(0);

    glDisable(GL_BLEND);
}

void Icon::render()
{
    for (int i = 0; i < positions.size(); i++)
    {
        update(i);
        renderPosition();
    }
}

void Icon::update(int positionIndex)
{
    // move the icons to their respective nodes
    if (positions.size() < 1)
        for (int i = 0; i < 4 /*size of arrays*/; i++)
            verts[i] = defaultVerts[i];
    else
        for (int i = 0; i < 4 /*size of arrays*/; i++)
            verts[i] = defaultVerts[i] + positions[positionIndex];

    // each icon only has 2 triangles and 4 uv's. 
    // the uv's never change so dont update them
    // the number of triangles doesnt change so just overwrite the buffer
    glBindVertexArray(vertexArray);

	glBindBuffer(GL_ARRAY_BUFFER, vertexBuffer);
	glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(verts[0]) * sizeof(verts), verts);
	glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 0, NULL);
	glEnableVertexAttribArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, 0);

	glBindVertexArray(0);
}

void Icon::loadImages()
{
	generateBuffer();

	image.clear();
	glDeleteTextures(1, &textureID);

	unsigned int error = lodepng::decode(image, imageWidth, imageHeight, filename);

	if (error)
		std::cout << "error " << error << ":" << lodepng_error_text(error) << std::endl;
	//creating 2D texture
	textureID = texture.generateTexture(image, imageWidth, imageHeight);
}