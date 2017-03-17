#include "Header.h"

// normalizing
void myNormalize(MyRay& ray)
{
	float len = length(ray.direction);
	ray.direction /= len;
}

void myNormalize(vec3 *ray)
{
	float len = length(*ray);
	*ray /= len;
}

float magnitude(vec3 ray)
{
	return sqrt((ray.x * ray.x) + (ray.y * ray.y + (ray.z + ray.z)));
}