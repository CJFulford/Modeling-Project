#include "Header.h"
#include <typeinfo>
#include <string>

using namespace glm;

// assuming that the incoming normal is normalized
vec3 Blinn_Phong(Ray *ray, float scalar, Object *object)
{
    #define SELECTED_BRIGHTEN_AMOUNT 1.5f
    vec3 intersect = ray->applyScalar(scalar),
        viewRay = normalize(ray->origin - intersect),
        lightRay = normalize(LIGHT_POS - intersect),
        halfRay = normalize((viewRay + lightRay)),
        normal = object->getNormal(intersect),
        col = object->colour * ((object->selected) ? SELECTED_BRIGHTEN_AMOUNT : 1.f);

    // specular colour is always white, add ambient for no black sides
	return  (col * AMBIENT) + 
            (WHITE * ((col * max(0.f, dot(normal, lightRay))) + 
                      (WHITE * pow(max(0.f, dot(normal, halfRay)), object->phong))));
}

vec3 getColour(Ray *ray, std::vector<Object*> *objectVec)
{
	for (unsigned int i = 0; i < objectVec->size(); i++)
	{
		Object *object = (*objectVec)[i];
		object->getVolume(ray);
	}

	float scalar = 0;
	Object *obj= ray->getClosestScalar(&scalar);

	if (scalar != 0)
		return Blinn_Phong(ray, scalar, obj);
	else 
		return BLUE * .1f; // navy blue background
}