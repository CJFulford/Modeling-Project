#include "Header.h"

vec3 shading(
	vec3 colourIn, 
	vec3 intersection, 
	vec3 origin, 
	vector<Light>& lightVec, 
	vec3 norm, 
	float phong, 
	vec3 specular)
{
	vec3 n = normalize(norm);
	vec3 colVec, ambient;

	vec3 v = normalize(origin - intersection);

	for (Light light : lightVec) {
		vec3 l = normalize(light.point - intersection);
		vec3 h = normalize((v + l));

		float diffuseFactor = std::max(0.f, dot(n, l));
		float specularFactor = std::max(0.f, dot(n, h));

		colVec += (colourIn * light.colour * diffuseFactor) + (specular * light.colour * pow(specularFactor, phong));
		ambient += light.ambient;
	}
	return colVec + (colourIn * ambient);
}

vec3 getColour(
	Ray& ray,
	vector<Sphere>& sphereVec,
	vector<Triangle>& triangleVec,
	vector<Light>& lightVec,
	int recursive)
{
	float scalar, sScalar, tScalar;
	vec3 colourVec, sCol, tCol;
	vec3 normal, sNorm, tNorm;
	vec3 specular, sSpecular, tSpecular;
	float phong, sPhong, tPhong;
	float reflect, sReflect, tReflect;

	sScalar = getNearestSphereScalar(ray, sphereVec, &sCol, &sNorm, &sPhong, &sSpecular, &sReflect);
	tScalar = getNearestTriangleScalar(ray, triangleVec, &tCol, &tNorm, &tPhong, &tSpecular, &tReflect);


	if (sScalar > 0 && (sScalar < tScalar || tScalar == 0))
	{
		normal = sNorm;
		colourVec = sCol;
		scalar = sScalar;
		phong = sPhong;
		specular = sSpecular;
		reflect = sReflect;
	}
	else if (tScalar > 0 && (tScalar < sScalar || sScalar == 0))
	{
		normal = tNorm;
		colourVec = tCol;
		scalar = tScalar;
		phong = tPhong;
		specular = tSpecular;
		reflect = tReflect;
	}
	else
		return BLACK;


	vec3 intersection = ray.origin + (scalar * ray.direction);

	if (recursive <= 0)
		return shading(colourVec, intersection, ray.origin, lightVec, normal, phong, specular);

	recursive--;
	Ray reflectedRay;
	reflectedRay.origin = intersection;
	reflectedRay.direction = normalize(ray.direction - (2.f * (dot(ray.direction, normal) * normal)));

	vec3 reflectedColourVec = getColour(reflectedRay, sphereVec, triangleVec, lightVec, recursive);

	return shading(colourVec + (reflect * reflectedColourVec), intersection, ray.origin, lightVec, normal, phong, specular);
}