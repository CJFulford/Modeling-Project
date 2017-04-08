#pragma once
#include <algorithm>
#include <glm\glm.hpp>
#include <vector>
#include <utility>
#include <glm/gtx/rotate_vector.hpp>
#include <complex>

// Mathematical values
#define PI				3.14159265359f
#define identity		glm::mat4(1.f)
#define FLOAT_ERROR 1e-4
#define radToDeg    (PI / 180)

#define SCALE_CHANGE .01f
#define MIN_SCALE .1f

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
    virtual void scale(bool enlarge) = 0;
    virtual void move(glm::vec3 move) = 0;
    virtual void rotate(glm::vec3 rotate) = 0;
    glm::vec3 center, colour;
    float phong, radius;
    bool selected = false, selectable = true;
};

// the actual objects
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
    void scale(bool enlarge) {}
    void move(glm::vec3 move) {}
    void rotate(glm::vec3 rotate) {}
    void select() {}
    void deselect() {}
    void breakBoolean(std::vector<Object*> *objectVec, int index) {}
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

struct Torus : Object
{
    Torus(glm::vec3 cent, float R, float r, glm::vec3 col, float ph) :
        R(R), r(r)
    {
        center = cent;
        radius = r;
        colour = col;
        phong = ph;
    }
    float R, r;


    const double TOLERANCE = 1.0e-8;
    // Returns n=0..numComplexValues, and fills in outArray with the n values from
    // inArray that are real-valued (i.e., whose imaginary parts are within TOLERANCE of 0.)
    // outArray must be large enough to receive numComplexValues values.
    int FilterRealNumbers(int numComplexValues, const std::complex<double> inArray[], double outArray[])
    {
        int numRealValues = 0;
        for (int i = 0; i < numComplexValues; ++i)
        {
            if (fabs(inArray[i].imag()) < TOLERANCE)
            {
                outArray[numRealValues++] = inArray[i].real();
            }
        }
        return numRealValues;
    }
    bool IsZero(std::complex<double> x)
    {
        return (fabs(x.real()) < TOLERANCE) && (fabs(x.imag()) < TOLERANCE);
    }
    // Returns n=0..2, the number of distinct real roots found for the equation
    //
    //     ax^2 + bx + c = 0
    //
    // Stores the roots in the first n slots of the array 'roots'.
    int SolveQuadraticEquation(std::complex<double> a, std::complex<double> b, std::complex<double> c, std::complex<double> roots[2])
    {
        if (IsZero(a))
        {
            if (IsZero(b))
            {
                // The equation devolves to: c = 0, where the variable x has vanished!
                return 0;   // cannot divide by zero, so there is no solution.
            }
            else
            {
                // Simple linear equation: bx + c = 0, so x = -c/b.
                roots[0] = -c / b;
                return 1;   // there is a single solution.
            }
        }
        else
        {
            const std::complex<double> radicand = b*b - 4.0*a*c;
            if (IsZero(radicand))
            {
                // Both roots have the same value: -b / 2a.
                roots[0] = -b / (2.0 * a);
                return 1;
            }
            else
            {
                // There are two distinct real roots.
                const std::complex<double> r = sqrt(radicand);
                const std::complex<double> d = 2.0 * a;

                roots[0] = (-b + r) / d;
                roots[1] = (-b - r) / d;
                return 2;
            }
        }
    }

    std::complex<double> cbrt(std::complex<double> a, int n)
    {
        /*
        This function returns one of the 3 complex cube roots of the complex number 'a'.
        The value of n=0..2 selects which root is returned.
        */

        const double TWOPI = 2.0 * 3.141592653589793238462643383279502884;

        double rho = pow(abs(a), 1.0 / 3.0);
        double theta = ((TWOPI * n) + arg(a)) / 3.0;
        return std::complex<double>(rho * cos(theta), rho * sin(theta));
    }

    // Returns n=0..3, the number of distinct real roots found for the equation
    //
    //     ax^3 + bx^2 + cx + d = 0
    //
    // Stores the roots in the first n slots of the array 'roots'.
    int SolveCubicEquation(std::complex<double> a, std::complex<double> b, std::complex<double> c, std::complex<double> d, std::complex<double> roots[3])
    {
        if (IsZero(a))
        {
            return SolveQuadraticEquation(b, c, d, roots);
        }

        b /= a;
        c /= a;
        d /= a;

        std::complex<double> S = b / 3.0;
        std::complex<double> D = c / 3.0 - S*S;
        std::complex<double> E = S*S*S + (d - S*c) / 2.0;
        std::complex<double> Froot = sqrt(E*E + D*D*D);
        std::complex<double> F = -Froot - E;

        if (IsZero(F))
        {
            F = Froot - E;
        }

        for (int i = 0; i < 3; ++i)
        {
            const std::complex<double> G = cbrt(F, i);
            roots[i] = G - D / G - S;
        }

        return 3;
    }

    // Returns n=0..4, the number of distinct real roots found for the equation
    //
    //     ax^4 + bx^3 + cx^2 + dx + e = 0
    //
    // Stores the roots in the first n slots of the array 'roots'.
    int SolveQuarticEquation(std::complex<double> a, std::complex<double> b, std::complex<double> c, std::complex<double> d, std::complex<double> e, std::complex<double> roots[4])
    {
        if (IsZero(a))
        {
            return SolveCubicEquation(b, c, d, e, roots);
        }

        // See "Summary of Ferrari's Method" in http://en.wikipedia.org/wiki/Quartic_function

        // Without loss of generality, we can divide through by 'a'.
        // Anywhere 'a' appears in the equations, we can assume a = 1.
        b /= a;
        c /= a;
        d /= a;
        e /= a;

        std::complex<double> b2 = b * b;
        std::complex<double> b3 = b * b2;
        std::complex<double> b4 = b2 * b2;

        std::complex<double> alpha = (-3.0 / 8.0)*b2 + c;
        std::complex<double> beta = b3 / 8.0 - b*c / 2.0 + d;
        std::complex<double> gamma = (-3.0 / 256.0)*b4 + b2*c / 16.0 - b*d / 4.0 + e;

        std::complex<double> alpha2 = alpha * alpha;
        std::complex<double> t = -b / 4.0;

        if (IsZero(beta))
        {
            std::complex<double> rad = sqrt(alpha2 - 4.0*gamma);
            std::complex<double> r1 = sqrt((-alpha + rad) / 2.0);
            std::complex<double> r2 = sqrt((-alpha - rad) / 2.0);

            roots[0] = t + r1;
            roots[1] = t - r1;
            roots[2] = t + r2;
            roots[3] = t - r2;
        }
        else
        {
            std::complex<double> alpha3 = alpha * alpha2;
            std::complex<double> P = -(alpha2 / 12.0 + gamma);
            std::complex<double> Q = -alpha3 / 108.0 + alpha*gamma / 3.0 - beta*beta / 8.0;
            std::complex<double> R = -Q / 2.0 + sqrt(Q*Q / 4.0 + P*P*P / 27.0);
            std::complex<double> U = cbrt(R, 0);
            std::complex<double> y = (-5.0 / 6.0)*alpha + U;
            if (IsZero(U))
            {
                y -= cbrt(Q, 0);
            }
            else
            {
                y -= P / (3.0 * U);
            }
            std::complex<double> W = sqrt(alpha + 2.0*y);

            std::complex<double> r1 = sqrt(-(3.0*alpha + 2.0*y + 2.0*beta / W));
            std::complex<double> r2 = sqrt(-(3.0*alpha + 2.0*y - 2.0*beta / W));

            roots[0] = t + (W - r1) / 2.0;
            roots[1] = t + (W + r1) / 2.0;
            roots[2] = t + (-W - r2) / 2.0;
            roots[3] = t + (-W + r2) / 2.0;
        }

        return 4;
    }

    // The above solvers are generalized for complex-valued 
    // coefficients and roots.
    // Below are helper methods that accept real-valued 
    // coefficients and return the subset of roots that are real-valued.

    inline int SolveQuadraticEquation(double a, double b, double c, double roots[2])
    {
        std::complex<double> croots[2];

        const int numComplexRoots = SolveQuadraticEquation(
            a,
            b,
            c,
            croots);

        return FilterRealNumbers(numComplexRoots, croots, roots);
    }

    inline int SolveCubicEquation(double a, double b, double c, double d, double roots[3])
    {
        std::complex<double> croots[3];

        const int numComplexRoots = SolveCubicEquation(
            a,
            b,
            c,
            d,
            croots);

        return FilterRealNumbers(numComplexRoots, croots, roots);
    }

    inline int SolveQuarticEquation(double a, double b, double c, double d, double e, double roots[4])
    {
        std::complex<double> croots[4];

        const int numComplexRoots = SolveQuarticEquation(
            a,
            b,
            c,
            d,
            e,
            croots);

        return FilterRealNumbers(numComplexRoots, croots, roots);
    }



    /*
    Finds where the ray interset the trus and adds all volumes to the ray's volume array
    resources for intersection
    http://www.cosinekitty.com/raytrace/chapter13_torus.html

    quartic equation solver has been copied then modified from the code provided on the same website
    */
    void getVolume(Ray *ray)
    {
        double T = 4.f * R * R;
        double G = T * ((ray->direction.x * ray->direction.x) + (ray->direction.y * ray->direction.y));
        double H = 2.f * T * ((ray->origin.x * ray->direction.x) + (ray->origin.y * ray->direction.y));
        double I = T * ((ray->origin.x * ray->origin.x) + (ray->origin.y * ray->origin.y));
        double J = (ray->direction.x * ray->direction.x) + (ray->direction.y * ray->direction.y) + (ray->direction.z * ray->direction.z);
        double K = 2.f * dot(ray->origin, ray->direction);
        double L = (ray->origin.x * ray->origin.x) + (ray->origin.y * ray->origin.y) + (ray->origin.z * ray->origin.z) + (R * R) - (r * r);

        std::vector<float> scalars;
        double roots[4];
        const int numRealRoots = SolveQuarticEquation(
            J*J,                    // coefficient of u^4
            2.0*J*K,                // coefficient of u^3
            2.0*J*L + K*K - G,      // coefficient of u^2
            2.0*K*L - H,            // coefficient of u^1 = u
            L*L - I,                // coefficient of u^0 = constant term
            roots                  // receives 0..4 real solutions
        );

        // We need to keep only the real roots.
        // There can be significant roundoff error in quartic solver, 
        // so we have to tolerate more slop than usual.
        const double SURFACE_TOLERANCE = 1.0e-4;
        int numPositiveRoots = 0;
        for (int i = 0; i < numRealRoots; ++i)
        {
            // Compact the array...
            if (roots[i] > SURFACE_TOLERANCE)
            {
               scalars.push_back(roots[i]);
            }
        }




        std::sort(scalars.begin(), scalars.end());
        switch (scalars.size())
        {
        case(1):
            //ray->volumes.push_back(Volume(scalars[0], scalars[0], this));
            break;
        case(2):
            ray->volumes.push_back(Volume(scalars[0], scalars[1], this));
            break;
        case(3):
            //ray->volumes.push_back(Volume(scalars[0], scalars[2], this));
            break;
        case(4):
            ray->volumes.push_back(Volume(scalars[0], scalars[1], this));
            ray->volumes.push_back(Volume(scalars[2], scalars[3], this));
            break;
        default:
            break;
        }
    }
    glm::vec3 getNormal(glm::vec3 intersection)
    {
        glm::vec3 point = intersection - center;
        float a = 1.f - (R / sqrt(point.x*point.x + point.y*point.y));
        return normalize(glm::vec3(a*point.x, a*point.y, point.z));
    }
    void scale(bool enlarge) {}
    void move(glm::vec3 move) {}
    void rotate(glm::vec3 rotate) {}
    void select() { selected = true; }
    void deselect() { selected = false; }
    void breakBoolean(std::vector<Object*> *objectVec, int index) {}
};

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

        bool obj1Torus = false; // need to know order of objects if one i a torus
        if (tempRay.volumes.size() == 2) obj1Torus = true;

        object2->getVolume(&tempRay);

        // if ray does not collide with both in the intersection
        if (tempRay.volumes.size() < 2)
            return;

        float entrance = 0, exit = 0;
        // special case for torus
        if (tempRay.volumes.size() == 3)
        {
            // obj1 == torus
            if (obj1Torus)
            {
                if (tempRay.volumes[0].exit <= tempRay.volumes[2].entrance ||
                    tempRay.volumes[2].exit <= tempRay.volumes[0].entrance ||
                    tempRay.volumes[1].exit <= tempRay.volumes[2].entrance ||
                    tempRay.volumes[2].exit <= tempRay.volumes[1].entrance)
                    return;

                entrance = glm::max(tempRay.volumes[0].entrance, tempRay.volumes[2].entrance);
                exit = glm::min(tempRay.volumes[0].exit, tempRay.volumes[2].exit);
                ray->volumes.push_back(Volume(entrance, exit, this));

                entrance = glm::max(tempRay.volumes[1].entrance, tempRay.volumes[2].entrance);
                exit = glm::min(tempRay.volumes[1].exit, tempRay.volumes[2].exit);
                ray->volumes.push_back(Volume(entrance, exit, this));
            }
            else
            {
                if (tempRay.volumes[0].exit <= tempRay.volumes[1].entrance ||
                    tempRay.volumes[1].exit <= tempRay.volumes[0].entrance ||
                    tempRay.volumes[0].exit <= tempRay.volumes[2].entrance ||
                    tempRay.volumes[2].exit <= tempRay.volumes[0].entrance)
                    return;

                entrance = glm::max(tempRay.volumes[0].entrance, tempRay.volumes[1].entrance);
                exit = glm::min(tempRay.volumes[0].exit, tempRay.volumes[1].exit);
                ray->volumes.push_back(Volume(entrance, exit, this));

                entrance = glm::max(tempRay.volumes[0].entrance, tempRay.volumes[2].entrance);
                exit = glm::min(tempRay.volumes[0].exit, tempRay.volumes[2].exit);
                ray->volumes.push_back(Volume(entrance, exit, this));
            }
        }
        else 
        {
            if (tempRay.volumes[0].exit <= tempRay.volumes[1].entrance ||
                tempRay.volumes[1].exit <= tempRay.volumes[0].entrance)
                return;
            entrance = glm::max(tempRay.volumes[0].entrance, tempRay.volumes[1].entrance);
            exit = glm::min(tempRay.volumes[0].exit, tempRay.volumes[1].exit);
            ray->volumes.push_back(Volume(entrance, exit, this));
        }
    }
    glm::vec3 getNormal(glm::vec3 intersection)
    {
        if (distance(intersection, object1->center) - object1->radius <= FLOAT_ERROR)
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
            std::vector<float> entrances;
            std::vector<float> exits;
            for (Volume volume : tempRay.volumes)
            {
                entrances.push_back(volume.entrance);
                exits.push_back(volume.exit);
            }
            std::sort(entrances.begin(), entrances.end());
            std::sort(exits.begin(), exits.end());

            ray->volumes.push_back(Volume(entrances.front(), exits.back(), this));
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

        // to toggle special conditions for torus
        bool obj1Torus = false;
        if (tempRay.volumes.size() == 2) obj1Torus = true;


        object2->getVolume(&tempRay);


        // if there is only 1 volume here, then the ray only collided with A in A - B.
        // lazy evaluation so if the size is < 2, it will not check the other conditions, 
        // avoiding a seg fault
        if (tempRay.volumes.size() < 2 ||
            tempRay.volumes[0].exit <= tempRay.volumes[1].entrance ||
            tempRay.volumes[1].exit <= tempRay.volumes[0].entrance)
            ray->volumes.push_back(Volume(tempRay.volumes[0].entrance, tempRay.volumes[0].exit, this));
        else if (tempRay.volumes.size() == 2)
        {
            Volume vol0 = tempRay.volumes[0], vol1 = tempRay.volumes[1];
            if (vol0.exit < vol1.entrance ||
                vol1.exit < vol0.entrance)
                return;

            if (vol0.entrance < vol1.entrance)
                ray->volumes.push_back(Volume(vol0.entrance, vol1.entrance, this));
            else if (vol1.exit < vol0.exit)
                ray->volumes.push_back(Volume(vol1.exit, vol0.exit, this));
        }
        else if (tempRay.volumes.size() == 3)
        {
            if (obj1Torus)
            {
                if (tempRay.volumes[0].exit <= tempRay.volumes[2].entrance ||
                    tempRay.volumes[2].exit <= tempRay.volumes[0].entrance ||
                    tempRay.volumes[1].exit <= tempRay.volumes[2].entrance ||
                    tempRay.volumes[2].exit <= tempRay.volumes[1].entrance)
                    return;

                Volume  vol0 = tempRay.volumes[0], 
                        vol1 = tempRay.volumes[1],
                        vol2 = tempRay.volumes[2];

                // compare first torus volume intersection to second object intersection
                if (vol0.entrance < vol2.entrance)
                    ray->volumes.push_back(Volume(vol0.entrance, vol2.entrance, this));
                else if (vol2.exit <= vol0.exit)
                    ray->volumes.push_back(Volume(vol2.exit, vol0.exit, this));

                // compare second torus volume intersection to second object intersection
                if (vol1.entrance < vol2.entrance)
                    ray->volumes.push_back(Volume(vol1.entrance, vol2.entrance, this));
                else if (vol2.exit <= vol1.exit)
                    ray->volumes.push_back(Volume(vol2.exit, vol1.exit, this));
            }
            else
            {
                if (tempRay.volumes[0].exit <= tempRay.volumes[1].entrance ||
                    tempRay.volumes[1].exit <= tempRay.volumes[0].entrance ||
                    tempRay.volumes[0].exit <= tempRay.volumes[2].entrance ||
                    tempRay.volumes[2].exit <= tempRay.volumes[0].entrance)
                    return;

                Volume  vol0 = tempRay.volumes[0],
                        vol1 = tempRay.volumes[1],
                        vol2 = tempRay.volumes[2];

                // compare first object intersection to first torus volume intersection
                if (vol0.entrance < vol1.entrance)
                    ray->volumes.push_back(Volume(vol0.entrance, vol1.entrance, this));
                else if (vol1.exit <= vol0.exit)
                    ray->volumes.push_back(Volume(vol1.exit, vol0.exit, this));

                // compare second object intersection to second torus volume intersection
                if (vol0.entrance < vol2.entrance)
                    ray->volumes.push_back(Volume(vol0.entrance, vol2.entrance, this));
                else if (vol2.exit <= vol0.exit)
                    ray->volumes.push_back(Volume(vol2.exit, vol0.exit, this));
            }
        }
    }
    glm::vec3 getNormal(glm::vec3 intersection)
    {
        if (abs(distance(intersection, center) - object1->radius) <= FLOAT_ERROR)
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