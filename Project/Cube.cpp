#include "Structs.h"

struct Plane : Object
{
    Plane(glm::vec3 point, glm::vec3 cubeCent, glm::vec3 normal, glm::vec3 col, float ph) :
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
    void scale(bool enlarge) {}
    void move(glm::vec3 move)
    {
        center += move;
    }
    void rotate(glm::vec3 rotate)
    {
        center = glm::rotateZ(glm::rotateY(glm::rotateX(center, rotate.x), rotate.y), rotate.z);
        normal = normalize(glm::rotateZ(glm::rotateY(glm::rotateX(normal, rotate.x), rotate.y), rotate.z));
    }
    void select() { selected = true; }
    void deselect() { selected = false; }
    void breakBoolean(std::vector<Object*> *objectVec, int index) {}
};

struct Cube : Object
{
#define XAXIS glm::vec3(1.f, 0.f, 0.f)
#define YAXIS glm::vec3(0.f, 1.f, 0.f)
#define ZAXIS glm::vec3(0.f, 0.f, 1.f)

    // variables
    Plane* planes[6];
    glm::vec3 top, side, front;
    float radi;

    // constructors
    Cube(glm::vec3 c, float r, glm::vec3 col, float ph) : radi(r)
    {
        radius = r;
        center = c;
        top = YAXIS;
        side = XAXIS;
        front = ZAXIS;
        colour = col;
        phong = ph;

        planes[0] = new Plane(center + (r * side), c, side, col, ph);
        planes[1] = new Plane(center + (r * -side), c, -side, col, ph);
        planes[2] = new Plane(center + (r * top), c, top, col, ph);
        planes[3] = new Plane(center + (r * -top), c, -top, col, ph);
        planes[4] = new Plane(center + (r * front), c, front, col, ph);
        planes[5] = new Plane(center + (r * -front), c, -front, col, ph);
    }
    ~Cube()
    {
        for (Plane* plane : planes)
            delete plane;
    };

    void getVolume(Ray *ray)
    {
        std::vector<float> points;

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
            // if point is contained in the cube, add to point vec
            if (contained) points.push_back(scalar);
        }

        // ignore single point result
        if (points.size() < 1) return;

        // add intersections to ray
        std::sort(points.begin(), points.end());
        ray->volumes.push_back(Volume(points.front(), points.back(), this));
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
        return glm::vec3(0.f);
    }
    void scale(bool enlarge)
    {
        if (enlarge)
        {
            radius += SCALE_CHANGE;
            for (Plane *plane : planes)
                plane->center = center + (radius * plane->normal);
        }
        else
        {
            radius = glm::max(MIN_SCALE, radius - SCALE_CHANGE);
            for (Plane *plane : planes)
                plane->center = center + (radius * plane->normal);
        }
    }
    void move(glm::vec3 move)
    {
        center += move;
        for (Plane *plane : planes)
            plane->move(move);
    }
    void rotate(glm::vec3 rotate)
    {
        /*top = normalize(glm::rotateZ(glm::rotateY(glm::rotateX(YAXIS, rotation.x), rotation.y), rotation.z));
        side = normalize(rotateZ(glm::rotateY(glm::rotateX(XAXIS, rotation.x), rotation.y), rotation.z));
        front = normalize(cross(side, top));*/
        for (Plane *plane : planes)
            plane->rotate(rotate);
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