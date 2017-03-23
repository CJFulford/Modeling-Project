#include "Header.h"

// getting the scalars to shapes
float getNearestSphereScalar(
	Ray ray, 
	vector<Sphere>& sphereVec, 
	vec3 *colourVec, 
	vec3 *normal, 
	float *phong, 
	vec3 *specular, 
	float *reflect)
{
	vec3 colVec = BLACK, norm, spec;
	float rayScalar = 0, tempScalar = 0, tempPhong = 0, tempReflect = 0;
	for (Sphere sphere : sphereVec)
	{
		tempScalar = 0;
		float rootValue = pow(dot(ray.direction, (ray.origin - sphere.center)), 2.f)
						- pow(length(ray.origin - sphere.center), 2.f)
						+ pow(sphere.radius, 2.f);
		if (rootValue >= 0)
		{
			float first = -dot(ray.direction, (ray.origin - sphere.center));
			if (rootValue > 0)
			{
				float	pos = first + sqrt(rootValue),
						neg = first - sqrt(rootValue);
				if (abs(neg) < abs(pos))
					tempScalar = neg;
				else
					tempScalar = pos;
			}
			else
				tempScalar = first;
		}

		if ((rayScalar == 0 || abs(tempScalar) < abs(rayScalar)) && tempScalar > FLOAT_ERROR)
		{
			norm = ((ray.origin + (tempScalar * ray.direction)) - sphere.center) / sphere.radius;
			rayScalar = tempScalar;
			colVec = sphere.colour;
			tempPhong = sphere.phong;
			spec = sphere.specular;
			tempReflect = sphere.reflect;
		}
	}
	*normal = normalize(norm);
	*colourVec = colVec;
	*phong = tempPhong;
	*specular = spec;
	*reflect = tempReflect;
	return rayScalar;
}


float getNearestTriangleScalar(
	Ray ray, 
	vector<Triangle>& triangleVec, 
	vec3 *colourVec, 
	vec3 *normal, 
	float *phong, 
	vec3 *specular, 
	float *reflect)
{
	vec3 colVec = BLACK, norm, spec;
	float rayScalar = 0, tempScalar = 0, tempPhong = 0, tempReflect = 0;
	for (Triangle tri : triangleVec) {
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

		if (((rayScalar == 0) || abs(tempScalar) < abs(rayScalar)) && tempScalar > FLOAT_ERROR)
		{
			tempPhong = tri.phong;
			norm = cross(tri.p3 - tri.p2, tri.p1 - tri.p2);
			rayScalar = tempScalar;
			colVec = tri.colour;
			spec = tri.specular;
			tempReflect = tri.reflect;
		}
	}
	*normal = normalize(norm);
	*phong = tempPhong;
	*colourVec = colVec;
	*specular = spec;
	*reflect = tempReflect;
	return rayScalar;
}