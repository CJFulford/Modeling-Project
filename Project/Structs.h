#pragma once
#include <algorithm>
#include <glm\glm.hpp>
#include <vector>
#include <utility>
#include <glm/gtx/rotate_vector.hpp>

// Mathematical values
#define PI				3.14159265359f
#define identity		glm::mat4(1.f)
#define FLOAT_ERROR 1e-4
#define radToDeg    (PI / 180)

#define SCALE_CHANGE .01f
#define MIN_SCALE .1f

struct Object;


struct Volume
{
	Volume(float entrance, float exit, Object *object) :
        entrance(entrance), exit(exit), object(object) {}
	float entrance, exit;
	Object *object;
};
struct Ray
{
	Ray(glm::vec3 orig, glm::vec3 dir) :
        origin(orig)
    {
        direction = normalize(dir);
    }
	glm::vec3 origin, direction;
	std::vector<Volume> volumes;

	Object * getClosestScalar(float *scalar)
	{
		float min = 0;
		Object *object = NULL;
		for (unsigned int i = 0; i < volumes.size(); i++)
		{
			Volume *volume = &volumes[i];

			if (volume->entrance > 0 && (volume->entrance < min || min == 0))
			{
				min = volume->entrance;
				object = volume->object;
			}
			else if (volume->exit > 0 && (volume->exit < min || min == 0))
			{
				min = volume->exit;
				object = volume->object;
			}
		}
		*scalar = min;
		return object;
	}
    glm::vec3 applyScalar(float scalar)
    {
        return origin + (scalar * direction);
    }
};

// the actual objects
struct Object
{
    glm::vec3 center = glm::vec3(0.f);
    glm::vec3 colour = glm::vec3(1.f);
    float phong = 0.f;
    float radius = 1.f;
    bool selected = false;
    bool selectable = true;

    virtual void getVolume(Ray *ray) = 0;
    virtual glm::vec3 getNormal(glm::vec3 intersection) = 0;
    virtual void breakBoolean(std::vector<Object*> *objectVec, int index) = 0;
    virtual void select() = 0;
    virtual void deselect() = 0;
    virtual void scale(bool enlarge) = 0;
    virtual void move(glm::vec3 move) = 0;
    virtual void rotate(glm::vec3 rotate) = 0;
};
struct Triangle;
struct Sphere;
struct Cube;
struct Torus;

// the combinations
struct Intersection;
struct Union;
struct Difference;


