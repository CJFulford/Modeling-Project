#include "Structs.h"

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
        if (object1->getNormal(intersection) != glm::vec3(0.f))
            return object1->getNormal(intersection);
        else
            return -object2->getNormal(intersection);
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