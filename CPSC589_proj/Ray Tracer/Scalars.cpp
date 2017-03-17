#include "Header.h"

// getting the scalars to shapes
float getNearestSphereScalar(
	MyRay ray, 
	vector<MySphere>& sphereVec, 
	vec3 *colourVec, 
	vec3 *normal, 
	float *phong, 
	vec3 *specular, 
	float *reflect)
{
	vec3 colVec = BLACK, norm, spec;
	float rayScalar = 0, tempScalar = 0, tempPhong = 0, tempReflect = 0;
	for (MySphere sphere : sphereVec)
	{
		tempScalar = 0;
		float rootValue = pow(dot(ray.direction, (ray.origin - sphere.center)), 2.f)
			- pow(length(ray.origin - sphere.center), 2.f)
			+ pow(sphere.radius, 2.f);
		if (rootValue < 0)
			tempScalar = 0;
		else
		{
			float first = -dot(ray.direction, (ray.origin - sphere.center));
			if (rootValue > 0)
			{
				float pos = first + sqrt(rootValue);
				float neg = first - sqrt(rootValue);
				if (abs(neg) < abs(pos))
					tempScalar = neg;
				else
					tempScalar = pos;
			}
			else
				tempScalar = first;
		}

		if (((rayScalar == 0) || abs(tempScalar) < abs(rayScalar)) && tempScalar > error)
		{
			norm = ((ray.origin + (tempScalar * ray.direction)) - sphere.center) / sphere.radius;
			rayScalar = tempScalar;
			colVec = sphere.colour;
			tempPhong = sphere.phong;
			spec = sphere.specular;
			tempReflect = sphere.reflect;
		}
	}
	myNormalize(&norm);
	*normal = norm;
	*colourVec = colVec;
	*phong = tempPhong;
	*specular = spec;
	*reflect = tempReflect;
	return rayScalar;
}

float getNearestSphereScalar(
	MyRay ray, 
	vector<MySphere>& sphereVec)
{
	float temp;
	vec3 tempVec;
	return getNearestSphereScalar(ray, sphereVec, &tempVec, &tempVec, &temp, &tempVec, &temp);
}

float getNearestTriangleScalar(
	MyRay ray, 
	vector<MyTriangle>& triangleVec, 
	vec3 *colourVec, 
	vec3 *normal, 
	float *phong, 
	vec3 *specular, 
	float *reflect)
{
	vec3 colVec = BLACK, norm, spec;
	float rayScalar = 0, tempScalar = 0, tempPhong = 0, tempReflect = 0;
	for (MyTriangle tri : triangleVec) {
		tempScalar = 0;
		float a = tri.a;
		float b = tri.b;
		float c = tri.c;
		float d = tri.d;
		float e = tri.e;
		float f = tri.f;
		float g = ray.direction.x;
		float h = ray.direction.y;
		float i = ray.direction.z;
		float j = tri.p1.x - ray.origin.x;
		float k = tri.p1.y - ray.origin.y;
		float l = tri.p1.z - ray.origin.z;

		float ei_hf = (e * i) - (h * f);
		float gf_di = (g * f) - (d * i);
		float dh_eg = (d * h) - (e * g);
		float ak_jb = (a * k) - (j * b);
		float jc_al = (j * c) - (a * l);
		float bl_kc = (b * l) - (k * c);

		float M = (a * ei_hf) + (b * gf_di) + (c * dh_eg);

		float tempScalar = -(((f * ak_jb) + (e * jc_al) + (d * bl_kc)) / M);

		float gamma = ((i * ak_jb) + (h * jc_al) + (g * bl_kc)) / M;
		if (gamma < 0 || gamma > 1) continue;

		float beta = ((j * ei_hf) + (k * gf_di) + (l * dh_eg)) / M;
		if (beta < 0 || beta  > 1 - gamma) continue;

		if (((rayScalar == 0) || abs(tempScalar) < abs(rayScalar)) && tempScalar > error)
		{
			tempPhong = tri.phong;
			norm = cross(tri.p3 - tri.p2, tri.p1 - tri.p2);
			rayScalar = tempScalar;
			colVec = tri.colour;
			spec = tri.specular;
			tempReflect = tri.reflect;
		}
	}
	myNormalize(&norm);
	*normal = norm;
	*phong = tempPhong;
	*colourVec = colVec;
	*specular = spec;
	*reflect = tempReflect;
	return rayScalar;
}

float getNearestTriangleScalar(
	MyRay ray, 
	vector<MyTriangle>& triangleVec)
{
	float temp;
	vec3 tempVec;
	return getNearestTriangleScalar(ray, triangleVec, &tempVec, &tempVec, &temp, &tempVec, &temp);
}

float getNearestPlaneScalar(
	MyRay ray, 
	vector<MyPlane>& planeVec, 
	vec3 *colourVec, 
	vec3 *normal, 
	float *phong, 
	vec3 *specular, 
	float *reflect)
{
	vec3 colVec = BLACK, norm, spec;
	float rayScalar = 0, tempScalar = 0, tempPhong = 0, tempReflect = 0;
	for (MyPlane plane : planeVec) {
		tempScalar = 0;

		float denominator = dot(ray.direction, plane.normal);
		float numerator = dot((plane.point - ray.origin), plane.normal);

		if (denominator == 0) continue;

		tempScalar = numerator / denominator;

		if (((rayScalar == 0) || abs(tempScalar) < abs(rayScalar)) && tempScalar > error)
		{
			norm = plane.normal;
			rayScalar = tempScalar;
			colVec = plane.colour;
			tempPhong = plane.phong;
			spec = plane.specular;
			tempReflect = plane.reflect;
		}
	}
	myNormalize(&norm);
	*normal = norm;
	*phong = tempPhong;
	*colourVec = colVec;
	*specular = spec;
	*reflect = tempReflect;
	return rayScalar;
}