#include "Structs.h"

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
    void scale(bool enlarge)
    {
        object1->scale(enlarge);
        object2->scale(enlarge);
    }
    void move(glm::vec3 move)
    {
        object1->move(move);
        object2->move(move);
    }
    void rotate(glm::vec3 rotate)
    {
        object1->rotate(rotate);
        object2->rotate(rotate);
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