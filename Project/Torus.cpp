#include "Header.h"

struct Torus : Object
{
    Torus(glm::vec3 cent, float R, float r, glm::vec3 col, float ph) :
        R(R), r(r)
    {
        center = cent;
        colour = col;
        phong = ph;
    }
    float R, r;

    /*
    Finds where the ray interset the trus and adds all volumes to the ray's volume array
    resources for intersection
    http://www.cosinekitty.com/raytrace/chapter13_torus.html
    */
    void getVolume(Ray *ray)
    {
        float G = 4.f * (R * R) * ((ray->origin.x * ray->origin.x) + (ray->origin.y * ray->origin.y));
        float H = 8.f * (R * R) * ((ray->direction.x * ray->origin.x) + (ray->direction.x * ray->origin.y));
        float I = 4.f * (R * R) * ((ray->direction.x * ray->direction.x) + (ray->direction.y * ray->direction.y));
        float J = (ray->origin.x + ray->origin.y + ray->origin.z) * (ray->origin.x + ray->origin.y + ray->origin.z);
        float K = 2.f * dot(ray->direction, ray->origin);
        float L = ((ray->direction.x + ray->direction.y + ray->direction.z) * (ray->direction.x + ray->direction.y + ray->direction.z)) + ((R * R) - (r * r));


    }
    glm::vec3 getNormal(glm::vec3 intersection)
    {
        return glm::vec3(0);
    }
    void scale(bool enlarge) {}
    void move(glm::vec3 move) {}
    void rotate(glm::vec3 rotate) {}
    void select() { selected = true; }
    void deselect() { selected = false; }
    void breakBoolean(std::vector<Object*> *objectVec, int index) {}
};