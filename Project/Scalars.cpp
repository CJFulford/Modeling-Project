#include "Header.h"

using namespace glm;

/*
Finds the scalar along the ray direction of the sphere closest to the ray origin.
Returns the scalar.
when function returns, normal and index will contain the normal of the sphere at
the intersection point and the sphereVec index of the sphere respectivly.

Algorithm from "Fundamentals of Computer Graphics 3rd ed. - P. Shirley, S. Marschner (CRC, 2009), p. 76
*/
float getSphereScalar(Ray *ray, std::vector<Sphere> *sphereVec, vec3 *normal, unsigned int *index)
{
	vec3 norm;
	float scalar = 0, tempScalar = 0;
    
    for (unsigned int i = 0; i < sphereVec->size(); i++)
	{
        Sphere *sphere = &(*sphereVec)[i];

		tempScalar = 0;

        // algorithm derives to the quadratic equation
        // value under the square root in the quadratic equation
		float rootValue = pow(dot(ray->direction, (ray->origin - sphere->center)), 2.f)
						- pow(length(ray->origin - sphere->center), 2.f)
						+ (sphere->radius * sphere->radius);


		if (rootValue >= 0)
		{
			float first = -dot(ray->direction, (ray->origin - sphere->center));
			if (rootValue > 0)
			{
				float	pos = first + sqrt(rootValue),
						neg = first - sqrt(rootValue);
				(abs(neg) < abs(pos)) ? tempScalar = neg : tempScalar = pos;
			}
			else
				tempScalar = first;
		}

        /* 
        If it is the first case 
        OR
        (If the new scalar is absolutly closer to the camera 
        AND 
        the reflection is not intersecting with itself) 
        */
		if ((scalar == 0 || abs(tempScalar) < abs(scalar)) && tempScalar > FLOAT_ERROR)
		{
			norm = ((ray->origin + (tempScalar * ray->direction)) - sphere->center);
			scalar = tempScalar;
            *index = i;
		}
	}
	*normal = normalize(norm);
	return scalar;
}

/*
Finds the scalar along the ray direction of the triangle closest to the ray origin.
Returns the scalar.
When function returns, normal and index will contain the normal of the triangle 
and the triangleVec index of the sphere respectivly.

Algorithm from "Fundamentals of Computer Graphics 3rd ed. - P. Shirley, S. Marschner (CRC, 2009), p. 77"
*/
float getTriScalar(Ray *ray, std::vector<Triangle> *triangleVec, unsigned int *index)
{
    vec3 norm;
	float scalar = 0, tempScalar = 0;

    // algorithm is Cramers rule that has been removed from matrix form for efficiency
    // this looks bad but it is a direct copy from the textbook and works
    for (unsigned int iter = 0; iter < triangleVec->size(); iter++) {
        Triangle *tri = &(*triangleVec)[iter];
		tempScalar = 0;
		float   a = tri->a,
		        b = tri->b,
		        c = tri->c,
                d = tri->d,
                e = tri->e,
                f = tri->f,
                g = ray->direction.x,
                h = ray->direction.y,
                i = ray->direction.z,
                j = tri->p1.x - ray->origin.x,
                k = tri->p1.y - ray->origin.y,
                l = tri->p1.z - ray->origin.z,

                ei_hf = (e * i) - (h * f),
                gf_di = (g * f) - (d * i),
                dh_eg = (d * h) - (e * g),
                ak_jb = (a * k) - (j * b),
                jc_al = (j * c) - (a * l),
                bl_kc = (b * l) - (k * c),

                M = (a * ei_hf) + (b * gf_di) + (c * dh_eg),

        // if gamma < 0 || gamma > 1 then the ray is missing the triangle
        gamma = ((i * ak_jb) + (h * jc_al) + (g * bl_kc)) / M;
		if (gamma < 0 || gamma > 1) continue;

        // if beta < 0 || beta > 1 then the ray is missing the triangle
		float beta = ((j * ei_hf) + (k * gf_di) + (l * dh_eg)) / M;
		if (beta < 0 || beta  > 1 - gamma) continue;
        
        // calculate the tempScalar if we need to
        tempScalar = -(((f * ak_jb) + (e * jc_al) + (d * bl_kc)) / M);

        /*
        If it is the first case
        OR
        (If the new scalar is absolutly closer to the camera
        AND
        the reflection is not intersecting with itself)
        */
		if (((scalar == 0) || abs(tempScalar) < abs(scalar)) && tempScalar > FLOAT_ERROR)
		{
			scalar = tempScalar;
            *index = iter;
		}
	}
	return scalar;
}