#pragma once
#include <glm\glm.hpp>
#include <vector>
#define FLOAT_ERROR 1e-6


struct Object;



struct Volume
{
	Volume(float entrance, float exit, Object *object) :
        entrance(entrance), exit(exit), object(object) {}
	float entrance, exit;
	Object *object;
};

struct Ray
{
	Ray(glm::vec3 orig, glm::vec3 dir) :
		origin(orig), direction(dir) {}
	glm::vec3 origin, direction;
	std::vector<Volume> volumes;

	Object * getClosestScalar(float *scalar)
	{
		float min = 0;
		Object *object = NULL;
		for (unsigned int i = 0; i < volumes.size(); i++)
		{
			Volume *volume = &volumes[i];

			if (volume->entrance > 0 && (volume->entrance < min || min == 0))
			{
				min = volume->entrance;
				object = volume->object;
			}
			else if (volume->exit > 0 && (volume->exit < min || min == 0))
			{
				min = volume->exit;
				object = volume->object;
			}
		}
		*scalar = min;
		return object;
	}
};



struct Object
{
	virtual void getVolume(Ray *ray) = 0;
	virtual glm::vec3 getNormal(glm::vec3 intersection) = 0;
	glm::vec3 colour;
	float phong;
};

struct Plane : Object
{
    Plane(glm::vec3 point, glm::vec3 normal) : point(point), normal(normalize(normal)) {}
    glm::vec3 point, normal;

    // returns scalar if it exists, returns null if ray and plane are parallel
    float getIntersection(Ray *ray)
    {
        float denominator = dot(ray->direction, normal);
        if (denominator == 0) return NULL;
        return dot((point - ray->origin), normal) / denominator;
    }

    void getVolume(Ray *ray) {}
    glm::vec3 getNormal(glm::vec3 intersection)
    {
        return normal;
    }
};

struct Triangle : Object
{
    Triangle(glm::vec3 p1, glm::vec3 p2, glm::vec3 p3, glm::vec3 col, float ph) :
        p1(p1), p2(p2), p3(p3)
    {
        a = p1.x - p2.x;
        b = p1.y - p2.y;
        c = p1.z - p2.z;

        d = p1.x - p3.x;
        e = p1.y - p3.y;
        f = p1.z - p3.z;

        normal = normalize(cross(p3 - p2, p1 - p2));
        colour = col;
        phong = ph;
    }
    glm::vec3 p1, p2, p3, normal;
    float a, b, c, d, e, f;

    /*
    Finds intersection point with ray and adds "volume" (infinitly thin so enterence and exit are the same) to the ray's volume array
    Algorithm from "Fundamentals of Computer Graphics 3rd ed. - P. Shirley, S. Marschner (CRC, 2009), p. 77"
    */
    void getVolume(Ray *ray)
    {
        // algorithm is Cramers rule that has been removed from matrix form for efficiency
        float	g = ray->direction.x,
            h = ray->direction.y,
            i = ray->direction.z,
            j = p1.x - ray->origin.x,
            k = p1.y - ray->origin.y,
            l = p1.z - ray->origin.z,

            ei_hf = (e * i) - (h * f),
            gf_di = (g * f) - (d * i),
            dh_eg = (d * h) - (e * g),
            ak_jb = (a * k) - (j * b),
            jc_al = (j * c) - (a * l),
            bl_kc = (b * l) - (k * c),

            M = (a * ei_hf) + (b * gf_di) + (c * dh_eg);

        // if gamma < 0 || gamma > 1 then the ray is missing the triangle
        float gamma = ((i * ak_jb) + (h * jc_al) + (g * bl_kc)) / M;
        if (gamma < 0 || gamma > 1) return;

        // if beta < 0 || beta > 1 then the ray is missing the triangle
        float beta = ((j * ei_hf) + (k * gf_di) + (l * dh_eg)) / M;
        if (beta < 0 || beta  > 1 - gamma) return;


        float scalar = -(((f * ak_jb) + (e * jc_al) + (d * bl_kc)) / M);
        ray->volumes.push_back(Volume(scalar, scalar, this));
    }
    glm::vec3 getNormal(glm::vec3 intersection)
    {
        return normal;
    }
};

struct Sphere : Object
{
	Sphere(glm::vec3 center, float radius, glm::vec3 col, float ph) : 
		center(center), radius(radius)
	{
		colour = col;
		phong = ph;
	}
	glm::vec3 center;
	float radius;

	/*
	Finds where the ray enters and exits this sphere, if at all, and adds the points to the rays volume array
	Algorithm from https://www.cs.princeton.edu/courses/archive/fall00/cs426/lectures/raycast/sld013.htm
	*/
	void getVolume(Ray *ray)
	{
		glm::vec3 L = center - ray->origin;
		float tca = dot(L, ray->direction);

		if (tca >= 0)
		{
			float d2 = dot(L, L) - (tca * tca),
				r2 = radius * radius;
			if (d2 <= r2)
			{
				float thc = sqrt(r2 - d2),
					t1 = tca - thc,
					t2 = tca + thc;
				if (abs(t1) < abs(t2))
					ray->volumes.push_back(Volume(t1, t2, this));
			}
		}
	}
	glm::vec3 getNormal(glm::vec3 intersection)
	{
		return glm::normalize(intersection - center);
	}
};

struct Cube : Object
{
    Cube(glm::vec3 center, glm::vec3 top, glm::vec3 side, float radius, glm::vec3 col, float ph) :
        center(center), top(normalize(top)), side(normalize(side)), radius(radius)
    {
        front = normalize(glm::cross(side, top));
        colour = col;
        phong = ph;

        planes[0] = new Plane(center + (radius * side), side);
        planes[1] = new Plane(center + (radius * -side), -side);
        planes[2] = new Plane(center + (radius * top), top);
        planes[3] = new Plane(center + (radius * -top), -top);
        planes[4] = new Plane(center + (radius * front), front);
        planes[5] = new Plane(center + (radius * -front), -front);
    }
    ~Cube() 
    {
        for (Plane* plane : planes)
            delete plane;
    };

    Plane* planes[6];
    glm::vec3 center, top, side, front;
    float radius;

    void getVolume(Ray *ray)
    {
        std::vector<float> points;

        // 6 is the number of planes in planes
        for (Plane* plane : planes)
        {
            float scalar = plane->getIntersection(ray);
            // if parallel, the scalar will be NULL
            if (!scalar) continue;

            // intersect position relative to the cube center
            glm::vec3 relativeIntersect = (ray->origin + (scalar * ray->direction)) - center,
                        corner = radius * (top + side + front);

            // if the intersection does not exceed the bounds of the normals
            if (abs(relativeIntersect.x) <= abs(corner.x) + FLOAT_ERROR &&
                abs(relativeIntersect.y) <= abs(corner.y) + FLOAT_ERROR &&
                abs(relativeIntersect.z) <= abs(corner.z) + FLOAT_ERROR)
                    points.push_back(scalar);
        }

        // add intersections to ray
        if (points.size() == 1) ray->volumes.push_back(Volume(points[0], points[0], this));
        else if (points.size() >= 2)
        {
            float min = points[0], max = points[0];
            for (float point : points)
            {
                if (point < min) min = point;
                else if (point > max) max = point;
            }
            ray->volumes.push_back(Volume(min, max, this));
        }
    }
    glm::vec3 getNormal(glm::vec3 intersection)
    {
        glm::vec3 relative = intersection - center;
        if (abs(relative.x) > abs(relative.y) && abs(relative.x) > abs(relative.z))
        {
            if (relative.x >= 0) return side;
            else return -side;
        }
        else if (abs(relative.y) > abs(relative.z))
        {
            if (relative.y >= 0) return top;
            else return -top;
        }
        else
        {
            if (relative.z >= 0) return front;
            else return -front;
        }
    }
};

struct Torus : Object
{
    Torus(glm::vec3 center, float mainRadius, float subRadius, glm::vec3 col, float ph) :
        center(center), mainRadius(mainRadius), subRadius(subRadius)
    {
        colour = col;
        phong = ph;
    }
    glm::vec3 center;
    float mainRadius, subRadius;

    /*
    Finds where the ray interset the trus and adds all volumes to the ray's volume array
    resources for intersection
    https://www.cl.cam.ac.uk/teaching/1999/AGraphHCI/SMAG/node2.html#SECTION00023400000000000000
    http://www.cosinekitty.com/raytrace/chapter13_torus.html
    http://www.wseas.org/multimedia/journals/computers/2013/025705-201.pdf

    */
    void getVolume(Ray *ray)
    {
    }
    glm::vec3 getNormal(glm::vec3 intersection)
    {
        return glm::vec3(0);
    }
};


