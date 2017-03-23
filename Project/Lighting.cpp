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
	vector<Light>& lightVec)
{
	float scalar, sScalar, tScalar;
	vec3	colourVec, sCol, tCol,
			normal, sNorm, tNorm,
			specular, sSpecular, tSpecular;
	float	phong, sPhong, tPhong,
			reflect, sReflect, tReflect;

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

	return shading(colourVec, ray.origin + (scalar * ray.direction), ray.origin, lightVec, normal, phong, specular);
}