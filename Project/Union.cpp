#include "Structs.h"

struct Union : Object
{
    Union(Object *obj1, Object *obj2) : object1(obj1), object2(obj2)
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
        if (object1->getNormal(intersection) != glm::vec3(0.f))
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