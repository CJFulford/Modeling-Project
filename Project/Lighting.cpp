#include "Header.h"

using namespace glm;

// assuming that the incoming normal is normalized
vec3 Blinn_Phong(Ray *ray, float scalar, vec3 colour, vec3 normal, float phong)
{
    vec3 intersect = ray->origin + (scalar * ray->direction),
        viewRay = normalize(ray->origin - intersect),
        lightRay = normalize(LIGHT_POS - intersect),
        halfRay = normalize((viewRay + lightRay));

    // specular colour is always white, add ambient for no black sides
	return  (colour * AMBIENT) + 
            (WHITE * ((colour * glm::max(0.f, dot(normal, lightRay))) + 
                      (WHITE * pow(glm::max(0.f, dot(normal, halfRay)), phong))));
}

vec3 getColour(Ray *ray, std::vector<Sphere> *sphereVec, std::vector<Triangle> *triangleVec)
{
	float sScalar, tScalar;
    vec3 sNorm;
    unsigned int sIndex, tIndex;

    // find the closest sphere and triangle by scalar
	sScalar = getSphereScalar(ray, sphereVec, &sNorm, &sIndex);
	tScalar = getTriangleScalar(ray, triangleVec, &tIndex);

    /*
    If the sphere is infront of the camera
    AND
    (the closest sphere is closer than the closest triangle
    OR
    the closest triangle is the reflecting object)
    */
	if (sScalar > 0 && (sScalar < tScalar || tScalar == 0))
	{
        Sphere *s = &(*sphereVec)[sIndex];
        return Blinn_Phong(ray, sScalar, s->colour, sNorm, s->phong);
	}
    /*
    If the triangle is infront of the camera
    AND
    (the closest triangle is closer than the closest sphere
    OR
    the closest sphere is the reflecting object)
    */
	else if (tScalar > 0 && (tScalar < sScalar || sScalar == 0))
	{
        Triangle *t = &(*triangleVec)[tIndex];
        return Blinn_Phong(ray, tScalar, t->colour, t->normal, t->phong);
	}
    // no intersection
	else
		return WHITE;
}