#include "Header.h"

// see if point is in shadow
bool checkShadow(
	vec3 intersection,
	vector<MySphere>& sphereVec,
	vector<MyPlane>& planeVec,
	vector<MyTriangle>& triangleVec,
	vector<MyLight>& lightVec)

{
	bool inShadow = false;
	MyRay ray;
	ray.origin = intersection;
	for (MyLight light : lightVec)
	{
		ray.direction = light.point - ray.origin;
		myNormalize(ray);
		float sScalar = getNearestSphereScalar(ray, sphereVec);
		float tScalar = getNearestTriangleScalar(ray, triangleVec);
		if (sScalar > error)
		{
			vec3 point = ray.origin + (sScalar * ray.direction);
			if (distance(ray.origin, light.point) > distance(ray.origin, point))
				inShadow = true;
		}
		if (tScalar > error)
		{
			vec3 point = ray.origin + (tScalar * ray.direction);
			if (distance(ray.origin, light.point) > distance(ray.origin, point))
				inShadow = true;
		}
	}
	return inShadow;
}

vec3 shading(
	vec3 colourIn, 
	vec3 intersection, 
	vec3 origin, 
	vector<MyLight>& lightVec, 
	vec3 n, 
	float phong, 
	vec3 specular)
{
	myNormalize(&n);
	vec3 colVec, ambient;

	vec3 v = origin - intersection;
	myNormalize(&v);

	for (MyLight light : lightVec) {
		vec3 l = light.point - intersection;
		myNormalize(&l);
		//vec3 h = (v + l) / magnitude(v + l);
		vec3 h = (v + l);
		myNormalize(&h);

		float diffuseFactor = std::max(0.f, dot(n, l));
		float specularFactor = std::max(0.f, dot(n, h));

		colVec += (colourIn * light.colour * diffuseFactor) + (specular * light.colour * pow(specularFactor, phong));
		ambient += light.ambient;
	}
	return colVec + (colourIn * ambient);
}

vec3 getColour(
	MyRay& ray,
	vector<MySphere>& sphereVec,
	vector<MyTriangle>& triangleVec,
	vector<MyPlane>& planeVec,
	vector<MyLight>& lightVec,
	int recursive)
{
	float scalar, sScalar, tScalar, pScalar;
	vec3 colourVec, sCol, tCol, pCol;
	vec3 normal, sNorm, tNorm, pNorm;
	vec3 specular, sSpecular, tSpecular, pSpecular;
	float phong, sPhong, tPhong, pPhong;
	float reflect, sReflect, tReflect, pReflect;

	sScalar = getNearestSphereScalar(ray, sphereVec, &sCol, &sNorm, &sPhong, &sSpecular, &sReflect);
	tScalar = getNearestTriangleScalar(ray, triangleVec, &tCol, &tNorm, &tPhong, &tSpecular, &tReflect);
	pScalar = getNearestPlaneScalar(ray, planeVec, &pCol, &pNorm, &pPhong, &pSpecular, &pReflect);


	if (sScalar > 0 && (sScalar < tScalar || tScalar == 0) && (sScalar < pScalar || pScalar == 0))
	{
		normal = sNorm;
		colourVec = sCol;
		scalar = sScalar;
		phong = sPhong;
		specular = sSpecular;
		reflect = sReflect;
	}
	else if (tScalar > 0 && (tScalar < sScalar || sScalar == 0) && (tScalar < pScalar || pScalar == 0))
	{
		normal = tNorm;
		colourVec = tCol;
		scalar = tScalar;
		phong = tPhong;
		specular = tSpecular;
		reflect = tReflect;
	}
	else if (pScalar > 0 && (pScalar < sScalar || sScalar == 0) && (pScalar < tScalar || tScalar == 0))
	{
		normal = pNorm;
		colourVec = pCol;
		scalar = pScalar;
		phong = pPhong;
		specular = pSpecular;
		reflect = pReflect;
	}
	else
		return BLACK;


	vec3 intersection = ray.origin + (scalar * ray.direction);

	bool inShadow = checkShadow(intersection, sphereVec, planeVec, triangleVec, lightVec);

	if (recursive <= 0 && inShadow)
		return colourVec * lightVec[0].ambient;
	else if (recursive <= 0 && !inShadow)
		return shading(colourVec, intersection, ray.origin, lightVec, normal, phong, specular);

	recursive--;
	MyRay reflectedRay;
	reflectedRay.origin = intersection;
	reflectedRay.direction = ray.direction - (2.f * (dot(ray.direction, normal) * normal));
	myNormalize(reflectedRay);

	vec3 reflectedColourVec = getColour(reflectedRay, sphereVec, triangleVec, planeVec, lightVec, recursive);

	if (inShadow)
		return (colourVec + (specular * reflectedColourVec)) * lightVec[0].ambient;
	else
		return shading(colourVec + (reflect * reflectedColourVec), intersection, ray.origin, lightVec, normal, phong, specular);
}