#include "Structs.h"

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
        if (abs(distance(intersection, center) - radius) < FLOAT_ERROR)
            return glm::normalize(intersection - center);
        return glm::vec3(0.f);
    }
    void scale(bool enlarge)
    {
        if (enlarge)
            radius += SCALE_CHANGE;
        else
            radius = glm::max(MIN_SCALE, radius - SCALE_CHANGE);
    }
    void move(glm::vec3 move)
    {
        center += move;
    }
    void rotate(glm::vec3 rotate) {}
    void select() { selected = true; }
    void deselect() { selected = false; }
    void breakBoolean(std::vector<Object*> *objectVec, int index) {}
};