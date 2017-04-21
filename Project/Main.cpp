#include "Header.h"
#include "CSGtree.h"
#include "TList.h"
#include "ImageBuffer.h"
#include <iostream>
#include <ctime>
#include <omp.h>
#include <glm/gtx/rotate_vector.hpp>
#include <stdio.h>

#define rayl	-1
#define rayr	1
#define rayt	1
#define rayb	-1

// definitions of variables used in ray tracing


glm::vec3 camOrigin = DEF_CAM_POS;
std::vector<Object*> objectVec;

// default scene consists of 2 spheres and a pyramid
void addFloor(std::vector<Object*> *objectVec)
{
    objectVec->push_back(new Triangle(
        glm::vec3(-2.f, -2.f, -2.f),        //p1
        glm::vec3(-2.f, -2.f, 2.f),         //p2
        glm::vec3( 2.f, -2.f, -2.f),        //p3
        glm::vec3(0.1f, 0.8f, 0.9f),        //colour
        50.f));                             //phong
    objectVec->push_back(new Triangle(
        glm::vec3( 2.f, -2.f, -2.f),
        glm::vec3(-2.f, -2.f,  2.f),
        glm::vec3( 2.f, -2.f,  2.f),
        glm::vec3(0.1f, 0.8f, 0.9f),
        50.f));
}

int main(int argc, char *argv[])
{
    GLFWwindow* window = generateWindow();

    // load glad
    gladLoadGL();
	
	ImageBuffer imageBuffer;
	imageBuffer.Initialize();

    addFloor(&objectVec);
	

	// variable initialization

	time_t startTime = 0;
	int frames = 0;
	float	w = -(rayr / (float)tan(FOV / 2));

    TList tlist = TList();
	CSGtree csg = CSGtree();
    Icon tListSymbols[] = {
        Icon("icons/A.png", glm::vec2(-0.95, -0.5)),
        Icon("icons/B.png", glm::vec2(-0.95, -0.6)),
        Icon("icons/Union.png", glm::vec2(-0.95,-0.7)),
        Icon("icons/Intersection.png", glm::vec2(-0.95, -0.8)),
        Icon("icons/Difference.png", glm::vec2(-0.95, -0.9)) }; 
	
    // line icons dont move. only update once
    for (Icon symbol : tListSymbols)
        symbol.update(0);

	while (!glfwWindowShouldClose(window))
	{
        startTime = time(NULL);

        // rotateX, rotateY, zoom. Do here for single calc
        glm::vec3 origin =  rotateY(
                                rotateX(camOrigin, rotate_x)
                                , rotate_y) 
                            * zoom;

		#pragma omp parallel for schedule(dynamic)
		for (int i = 0; i < RENDER_WINDOW_WIDTH; i++)
		{
			for (int j = 0; j < RENDER_WINDOW_HEIGHT; j++)
			{
    			glm::vec3 colourVec = BLACK;
				float u = rayl + ((rayr - rayl)*(i + .5f)) / RENDER_WINDOW_WIDTH;
				float v = rayb + ((rayt - rayb)*(j + .5f)) / RENDER_WINDOW_HEIGHT;

                // construct the ray
                // rotate the direction along the x axis then the y axis
				Ray ray(origin, 
                        glm::normalize(
                            rotateY(
                                rotateX(glm::vec3(u, v, w), rotate_x)
                                , rotate_y)));

				colourVec = getColour(&ray, &objectVec);
				imageBuffer.SetPixel(i, j + (WHOLE_HEIGHT - RENDER_WINDOW_HEIGHT), colourVec);
			}
		}
        #pragma omp barrier
        imageBuffer.Render();

        // compute and draw the tlists
        tlist.getLines(&tlistRay);
        tlist.render();

        for (Icon symbol : tListSymbols)
            symbol.render();
		
		if (selected1 != -1)
		{
            if (selected2 != -1)
                csg.constructInfo(objectVec[selected1], objectVec[selected2], 1/*base level of tree*/);
            else
			    csg.constructInfo(objectVec[selected1], 1/*base level of tree*/);
            csg.update();
		}
			
		
		csg.render();
        glfwSwapBuffers(window);

        // update and print frames per second
		frames++;
		if (difftime(time(NULL), startTime) >= 1)
		{
            printf("\rFPS: %i", frames);
			frames = 0;
		}

        glfwPollEvents();
	}

	// clean up allocated resources before exit

	glfwDestroyWindow(window);
	glfwTerminate();

    for (Object *object : objectVec)
        delete object;
	return 0;
}
