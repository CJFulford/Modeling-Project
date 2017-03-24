#include "Header.h"

using namespace glm;

/*
Finds the scalar along the ray direction of the sphere closest to the ray origin.
Returns the scalar.
when function returns, normal and index will contain the normal of the sphere at
the intersection point and the sphereVec index of the sphere respectivly.

Algorithm from https://www.cs.princeton.edu/courses/archive/fall00/cs426/lectures/raycast/sld013.htm
*/
float getSphereScalar(Ray *ray, std::vector<Sphere> *sphereVec, vec3 *normal, unsigned int *index)
{
	vec3 norm;
	float scalar = 0;
    
    for (unsigned int i = 0; i < sphereVec->size(); i++)
	{
        Sphere *sphere = &(*sphereVec)[i];
		float tempScalar = 0;


        vec3 L = sphere->center - ray->origin;
        float tca = dot(L, ray->direction);

        if (tca >= 0)
        {
            float d2 = dot(L, L) - (tca * tca),
                    r2 = sphere->radius * sphere->radius;
            if (d2 <= r2)
            {
                float thc = sqrt(r2 - d2),
                        t1 = tca - thc,
                        t2 = tca + thc;
                (abs(t1) < abs(t2)) ? tempScalar = t1 : tempScalar = t2;
            }
        }

        /* 
        If it is the first case 
        OR
        (If the new scalar is absolutly closer to the camera 
        AND 
        the reflection is not intersecting with itself) 
        */
		if ((scalar == 0 || abs(tempScalar) < abs(scalar)) && tempScalar > 0)
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
When function returns, index will contain the normal of the triangle

Algorithm from "Fundamentals of Computer Graphics 3rd ed. - P. Shirley, S. Marschner (CRC, 2009), p. 77"
*/
float getTriangleScalar(Ray *ray, std::vector<Triangle> *triangleVec, unsigned int *index)
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
		if (((scalar == 0) || abs(tempScalar) < abs(scalar)) && tempScalar > 0)
		{
			scalar = tempScalar;
            *index = iter;
		}
	}
	return scalar;
}

/*
Finds the scalar along the ray direction of the torus closest to the ray origin.
Returns the scalar.
When function returns, normal and index will contain the normal of the torus
and the torusVec index of the torus respectivly.

resources for intersection
https://www.cl.cam.ac.uk/teaching/1999/AGraphHCI/SMAG/node2.html#SECTION00023400000000000000
http://www.cosinekitty.com/raytrace/chapter13_torus.html
http://www.wseas.org/multimedia/journals/computers/2013/025705-201.pdf

*/
float getTorusScalar(Ray *ray, std::vector<Torus> *torusVec, vec3 *normal, unsigned int *index)
{
    vec3 norm;
    float scalar = 0;

    for (unsigned int i = 0; i < torusVec->size(); i++)
    {
        Torus *torus = &(*torusVec)[i];
        float tempScalar = 0;





        /*
        If it is the first case
        OR
        (If the new scalar is absolutly closer to the camera
        AND
        the reflection is not intersecting with itself)
        */
        if ((scalar == 0 || abs(tempScalar) < abs(scalar)) && tempScalar > 0)
        {
            // TODO norm = ((ray->origin + (tempScalar * ray->direction)) - torus->center);
            scalar = tempScalar;
            *index = i;
        }
    }
    *normal = normalize(norm);
    return scalar;
}