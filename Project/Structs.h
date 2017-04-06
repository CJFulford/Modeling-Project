#pragma once
#include <algorithm>
#include <glm\glm.hpp>
#include <vector>
#include <utility>
#include <glm/gtx/rotate_vector.hpp>

// Mathematical values
#define PI				3.14159265359f
#define identity		glm::mat4(1.f)
#define FLOAT_ERROR 1e-4
#define radToDeg    (PI / 180)

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
        origin(orig)
    {
        direction = normalize(dir);
    }
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
    glm::vec3 applyScalar(float scalar)
    {
        return origin + (scalar * direction);
    }
};



struct Object
{
	virtual void getVolume(Ray *ray) = 0;
    virtual glm::vec3 getNormal(glm::vec3 intersection) = 0;
    virtual void breakBoolean(std::vector<Object*> *objectVec, int index) = 0;
    virtual void select() = 0;
    virtual void deselect() = 0;
	glm::vec3 center, colour;
	float phong, radius;
    bool selected = false, selectable = true;
};

// the actual objects
struct Plane : Object
{
    Plane(glm::vec3 point, glm::vec3 normal, glm::vec3 col, float ph) : 
        normal(normalize(normal)) 
    {
        center = point;
        colour = col;
        phong = ph;
    }
    glm::vec3 normal;

    // returns scalar if it exists, returns null if ray and plane are parallel
    float getIntersection(Ray *ray)
    {
        float denominator = dot(ray->direction, normal);
        if (denominator == 0) return NULL;
        return dot((center - ray->origin), normal) / denominator;
    }

    void getVolume(Ray *ray) {}
    glm::vec3 getNormal(glm::vec3 intersection)
    {
        return normal;
    }
    void select() { selected = true; }
    void deselect() { selected = false; }
    void breakBoolean(std::vector<Object*> *objectVec, int index){}
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
        selectable = false;
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
    void select() {}
    void deselect() {}
    void breakBoolean(std::vector<Object*> *objectVec, int index){}
};

struct Sphere : Object
{
	Sphere(glm::vec3 cent, float radi, glm::vec3 col, float ph) : radi(radi)
	{
        radius = radi;
        center = cent;
		colour = col;
		phong = ph;
	}
    float radi;

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
    void select() { selected = true; }
    void deselect() { selected = false; }
    void breakBoolean(std::vector<Object*> *objectVec, int index){}
};

struct Cube : Object
{
#define XAXIS glm::vec3(1.f, 0.f, 0.f)
#define YAXIS glm::vec3(0.f, 1.f, 0.f)
#define ZAXIS glm::vec3(0.f, 0.f, 1.f)

    // variables
    Plane* planes[6];
    glm::vec3 center, top, side, front;
    float radius, xRot, yRot;

    // constructors
    Cube(glm::vec3 c, float x, float y, float r, glm::vec3 col, float ph) :
        center(c), radius(r)
    {
        xRot = x * radToDeg;
        yRot = y * radToDeg;
        top = normalize(glm::rotateY(glm::rotateX(YAXIS, xRot), yRot));
        side = normalize(glm::rotateY(glm::rotateX(XAXIS, xRot), yRot));
        front = normalize(cross(side, top));
        colour = col;
        phong = ph;

        planes[0] = new Plane(center + (radius * side), side, col, ph);
        planes[1] = new Plane(center + (radius * -side), -side, col, ph);
        planes[2] = new Plane(center + (radius * top), top, col, ph);
        planes[3] = new Plane(center + (radius * -top), -top, col, ph);
        planes[4] = new Plane(center + (radius * front), front, col, ph);
        planes[5] = new Plane(center + (radius * -front), -front, col, ph);
    }
    ~Cube() 
    {
        for (Plane* plane : planes)
            delete plane;
    };

    void getVolume(Ray *ray)
    {
        std::vector<std::pair<float, Plane*>> points;

        // 6 is the number of planes in planes
        for (Plane *plane : planes)
        {
            float scalar = plane->getIntersection(ray);
            // if parallel, the scalar will be NULL
            if (!scalar) continue;


            bool contained = true;
            for (Plane *planeCheck : planes)
            {
                glm::vec3 relativeIntersect = normalize(ray->applyScalar(scalar) - planeCheck->center);
                if (dot(relativeIntersect, planeCheck->normal) > FLOAT_ERROR)
                {
                    contained = false;
                    break;
                }
            }

            if (!contained) continue;

            points.push_back(std::make_pair(scalar, plane));
        }

        // add intersections to ray
        if (points.size() == 2)
        {
            float p0 =  points[0].first,    p1 = points[1].first;
            Plane *i0 = points[0].second,   *i1 = points[1].second;
            ray->volumes.push_back(Volume(glm::min(p0, p1), 
                                          glm::max(p0, p1), 
                                          this));
        }
        else if (points.size() >= 2)
        {
            std::pair<float, Plane*> minP = points[0], maxP = points[0];
            for (unsigned int i = 1; i < points.size(); i++)
            {
                if (points[i].first < minP.first)
                    minP = points[i];
                else if (points[i].first > maxP.first)
                    maxP = points[i];
            }
            ray->volumes.push_back(Volume(minP.first, 
                                          maxP.first, 
                                          this));
        }
    }
    glm::vec3 getNormal(glm::vec3 intersection)
    {
        for (Plane *plane : planes)
        {
            glm::vec3 relativeIntersect = normalize(intersection - plane->center);
            if (abs(dot(relativeIntersect, plane->normal)) <= FLOAT_ERROR)
            {
                return plane->normal;
            }
        }
        return glm::vec3(1.f);
    }
    void select()
    {
        selected = true;
        for (Plane *plane : planes)
            plane->selected = true;
    }
    void deselect()
    {
        selected = false;
        for (Plane *plane : planes)
            plane->selected = false;
    }
    void breakBoolean(std::vector<Object*> *objectVec, int index) {}
};

struct Torus : Object
{
    Torus(glm::vec3 cent, float mainRadius, float subRadius, glm::vec3 col, float ph) :
         mainRadius(mainRadius), subRadius(subRadius)
    {
        center = cent;
        colour = col;
        phong = ph;
    }
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
    void select() { selected = true; }
    void deselect() { selected = false; }
    void breakBoolean(std::vector<Object*> *objectVec, int index){}
};

// the combinations
struct Intersection : Object
{
    Intersection(Object *obj1, Object *obj2) : object1(obj1), object2(obj2) 
    {
        phong = (obj1->phong + obj2->phong) / 2.f;
        colour = (obj1->colour + obj2->colour) / 2.f;
        center = (obj1->center + obj2->center) / 2.f;
        selected = false;
    }
    Object *object1, *object2;
    
    void getVolume(Ray *ray)
    {
        Ray tempRay(ray->origin, ray->direction);
        object1->getVolume(&tempRay);
        object2->getVolume(&tempRay);

        // if ray does not collide with both in the intersection
        if (tempRay.volumes.size() != 2 ||
            tempRay.volumes[0].exit <= tempRay.volumes[1].entrance ||
            tempRay.volumes[1].exit <= tempRay.volumes[0].entrance) 
                return;

        float entrance = glm::max(tempRay.volumes[0].entrance, tempRay.volumes[1].entrance);
        float exit = glm::min(tempRay.volumes[0].exit, tempRay.volumes[1].exit);
        ray->volumes.push_back(Volume(entrance, exit, this));
    }
    glm::vec3 getNormal(glm::vec3 intersection)
    {
        if (distance(intersection, object1->center) >= distance(intersection, object2->center))
            return object1->getNormal(intersection);
        else
            return object2->getNormal(intersection);
    }
    void select() { selected = true; }
    void deselect() { selected = false; }
    void breakBoolean(std::vector<Object*> *objectVec, int index)
    {
        objectVec->erase(objectVec->begin() + index);
        objectVec->push_back(object1);
        objectVec->push_back(object2);
        delete this;
    }
};

struct Union : Object
{
    Union(Object *obj1, Object *obj2) : object1(obj1), object2(obj2)
    {
        phong  = (obj1->phong + obj2->phong) / 2.f;
        colour = (obj1->colour + obj2->colour) / 2.f;
        center = (obj1->center + obj2->center) / 2.f;
        selected = false;
    }
    Object *object1, *object2;

    void getVolume(Ray *ray)
    {
        Ray tempRay(ray->origin, ray->direction);
        object1->getVolume(&tempRay);
        object2->getVolume(&tempRay);

        // ray needs to intersect both to collide with the intersection
        if (tempRay.volumes.size() == 1)
            ray->volumes.push_back(Volume(tempRay.volumes[0].entrance, tempRay.volumes[0].exit, this));
        else if (tempRay.volumes.size() > 1)
        {
            float entrance = glm::min(tempRay.volumes[0].entrance, tempRay.volumes[1].entrance);
            float exit = glm::max(tempRay.volumes[0].exit, tempRay.volumes[1].exit);
            ray->volumes.push_back(Volume(entrance, exit, this));
        }
    }
    glm::vec3 getNormal(glm::vec3 intersection)
    {
        if (distance(intersection, object1->center) <= distance(intersection, object2->center))
            return object1->getNormal(intersection);
        else
            return object2->getNormal(intersection);
    }
    void select() { selected = true; }
    void deselect() { selected = false; }
    void breakBoolean(std::vector<Object*> *objectVec, int index)
    {
        objectVec->erase(objectVec->begin() + index);
        objectVec->push_back(object1);
        objectVec->push_back(object2);
        delete this;
    }
};

struct Difference : Object
{
    Difference(Object *obj1, Object *obj2) : object1(obj1), object2(obj2)
    {
        phong = (obj1->phong + obj2->phong) / 2.f;
        colour = (obj1->colour + obj2->colour) / 2.f;
        center = (obj1->center + obj2->center) / 2.f;
        selected = false;
    }
    Object *object1, *object2;

    void getVolume(Ray *ray)
    {
        Ray tempRay(ray->origin, ray->direction);
        


        object1->getVolume(&tempRay);
        // if there is no collision with A in A - B then there is no collision with the difference
        if (tempRay.volumes.size() == 0) return;


        object2->getVolume(&tempRay);


        // if there is only 1 volume here, then the ray only collided with A in A - B.
        // lazy evaluation so if the size is < 2, it will not check the other conditoins, 
        // avoiding a seg fault
        if (tempRay.volumes.size() < 2 ||
            tempRay.volumes[0].exit <= tempRay.volumes[1].entrance ||
            tempRay.volumes[1].exit <= tempRay.volumes[0].entrance)
            ray->volumes.push_back(Volume(tempRay.volumes[0].entrance, tempRay.volumes[0].exit, this));
        else if (tempRay.volumes.size() == 2)
        {

            if (tempRay.volumes[0].exit <= tempRay.volumes[1].entrance ||
                tempRay.volumes[1].exit <= tempRay.volumes[0].entrance)
                return;



            Volume vol1 = tempRay.volumes[0], vol2 = tempRay.volumes[1];
            float entrance, exit;
            if (vol1.entrance < vol2.entrance)
            {
                entrance = vol1.entrance;
                exit = vol2.entrance;
                ray->volumes.push_back(Volume(entrance, exit, this));
            }
            else if (vol2.exit <= vol1.exit)
            {
                entrance = vol2.exit;
                exit = vol1.exit;
                ray->volumes.push_back(Volume(entrance, exit, this));
            }


        }
    }
    glm::vec3 getNormal(glm::vec3 intersection)
    {
        if (abs(distance(intersection, object1->center) - object1->radius) <= FLOAT_ERROR)
            return object1->getNormal(intersection);
        else
            return -object2->getNormal(intersection);
    }
    void select() { selected = true; }
    void deselect() { selected = false; }
    void breakBoolean(std::vector<Object*> *objectVec, int index)
    {
        objectVec->erase(objectVec->begin() + index);
        objectVec->push_back(object1);
        objectVec->push_back(object2);
        delete this;
    }
};