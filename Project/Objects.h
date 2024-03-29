#pragma once
#include "Tools.h"
#include <algorithm>
#include <glm/glm.hpp>
#include <vector>
#include <glm/gtx/rotate_vector.hpp>
#include <complex>

// Mathematical values
#define PI			3.14159265359f
#define identity	glm::mat4(1.f)
#define FLOAT_ERROR 1e-4
#define radToDeg    (PI / 180)
#define ZERO_VECTOR glm::vec3(0.f)

#define BASE_CENTER ZERO_VECTOR
#define BASE_RADIUS .5f
#define BASE_PHONG  50
#define BASE_COLOUR glm::vec3(1.f, .4f, 0.f)

#define SCALE_CHANGE .01f
#define MIN_SCALE .1f

#define XAXIS glm::vec3(1.f, 0.f, 0.f)
#define YAXIS glm::vec3(0.f, 1.f, 0.f)
#define ZAXIS glm::vec3(0.f, 0.f, 1.f)

extern float trotate_x;
extern float trotate_y;
extern float rotate_x;
extern float rotate_y;

// need this defined for object as the virtual function getVolume takes in a Ray
struct Ray;

// need this here as the ray's pushVolume funciton takes in an object
struct Object
{
    // object variables
    glm::vec3 center = ZERO_VECTOR;
    glm::vec3 colour = ZERO_VECTOR;
    glm::vec3 rotation = ZERO_VECTOR;
    float phong = 50.f;
    float radius = 1.f;
    bool selected = false;
    bool selectable = true;
    bool differenceB = false;

    // object derivative functions
    virtual void getVolume(Ray *ray) = 0;
    virtual glm::vec3 getNormal(glm::vec3 intersection) = 0;
    virtual void breakBoolean(std::vector<Object*> *objectVec, int index) = 0;
    virtual void scale(bool enlarge) = 0;
    virtual void move(glm::vec3 move) = 0;
    virtual void rotate(glm::vec3 rotate) = 0;

    // object defined functions
    void select() { selected = true; }
    void deselect() { selected = false; }
};

// ray
struct Volume
{
    Volume(float entrance, float exit, glm::vec3 entranceNormal, glm::vec3 exitNormal, Object *object) :
        entrance(entrance), exit(exit), object(object), entranceNormal(entranceNormal), exitNormal(exitNormal) {}
    float entrance, exit;
    glm::vec3 entranceNormal, exitNormal;
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

    Object * getClosestScalar(float *scalar, glm::vec3 *normal)
    {
        float min = 0;
        glm::vec3 norm(0.f);
        Object *object = NULL;
        for (unsigned int i = 0; i < volumes.size(); i++)
        {
            Volume *volume = &volumes[i];

            if (volume->entrance > 0 && (volume->entrance < min || min == 0))
            {
                min = volume->entrance;
                norm = volume->entranceNormal;
                object = volume->object;
            }
            else if (volume->exit > 0 && (volume->exit < min || min == 0))
            {
                min = volume->exit;
                norm = volume->exitNormal;
                object = volume->object;
            }
        }
        *normal = norm;
        *scalar = min;
        return object;

    }
    glm::vec3 applyScalar(float scalar)
    {
        return origin + (scalar * direction);
    }
    std::vector<Volume> getVolume() { return volumes; }
    void pushVolume(float s1, float s2, Ray *ray, Object *obj)
    {
        glm::vec3 norm1 = obj->getNormal(ray->applyScalar(s1));
        glm::vec3 norm2 = obj->getNormal(ray->applyScalar(s2));
        volumes.push_back(Volume(s1, s2, norm1, norm2, obj));
    }
    void pushVolumeBoolean(float s1, float s2, glm::vec3 normal1, glm::vec3 normal2, Object *obj)
    {
        volumes.push_back(Volume(s1, s2, normal1, -normal2, obj));
    }
    Ray generateTempRay(glm::vec3 center, glm::vec3 rotation)
    {
        return Ray(
            glm::rotateX(glm::rotateY(glm::rotateZ(origin - center, -rotation.z), -rotation.y), -rotation.x),
            glm::rotateX(glm::rotateY(glm::rotateZ(direction, -rotation.z), -rotation.y), -rotation.x));
    }
};

// uncreateable objects
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
        ray->pushVolume(scalar, scalar, ray, this);
    }
    glm::vec3 getNormal(glm::vec3 intersection)
    {
        return normal;
    }
    void scale(bool enlarge) {}
    void move(glm::vec3 move) {}
    void rotate(glm::vec3 rotate) {}
    void breakBoolean(std::vector<Object*> *objectVec, int index) {}
};
struct Plane : Object
{
    glm::vec3 normal;
    Plane(glm::vec3 point, glm::vec3 normal, glm::vec3 col) :
        normal(normalize(normal))
    {
        center = point;
        colour = col;
        phong = BASE_PHONG;
    }

    // returns scalar if it exists, returns null if ray and plane are parallel
    float getIntersection(Ray *ray)
    {
        float denominator = dot(ray->direction, normal);
        if (denominator == 0) return NULL;
        return dot((center - ray->origin), normal) / denominator;
    }
    void getVolume(Ray *ray) {}
    glm::vec3 getNormal(glm::vec3 intersection) { return (differenceB) ? -normal : normal; }
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
    void breakBoolean(std::vector<Object*> *objectVec, int index) {}
};
struct RayCylinder : Object
{
#define BASE_LENGTH .5f
    glm::vec3 rotation = ZERO_VECTOR;

    Ray *orientation;

    RayCylinder(Ray *ray)
    {
        orientation = ray;
        radius = .01f;
        center = BASE_CENTER;
        colour = ZERO_VECTOR;
        phong = BASE_PHONG;
    }

    void getVolume(Ray *ray)
    {
        glm::vec3 origin = rotateX(
            rotateY(ray->origin, -(trotate_y))
            , -(trotate_x + (PI / 2.f)));

        Ray tempRay(origin,
            glm::normalize(
                rotateX(
                    rotateY(ray->direction, -(trotate_y))
                    , -(trotate_x + (PI / 2.f)))));

        // at most 1 plane intersection. Need to check the cylinder now
        float A = (tempRay.direction.x * tempRay.direction.x) + (tempRay.direction.z * tempRay.direction.z);
        float B = 2.f * ((tempRay.origin.x * tempRay.direction.x) + (tempRay.origin.z * tempRay.direction.z));
        float C = (tempRay.origin.x * tempRay.origin.x) + (tempRay.origin.z * tempRay.origin.z) - (radius * radius);

        // now we solve the quadratic equation
        float rootTerm = (B * B) - (4.f * A * C);
        if (rootTerm < 0) return;   // no solutions, misses the cylinder

        rootTerm = sqrt(rootTerm);

        float scalars[2] = 
        { 
            (-B - rootTerm) / (2.f * A),
            (-B + rootTerm) / (2.f * A) 
        };

        ray->pushVolume(glm::min(scalars[0], scalars[1]), glm::min(scalars[0], scalars[1]), ray, this);
    }
    glm::vec3 getNormal(glm::vec3 intersection)
    {
        glm::vec3 norm(0.f);
        glm::vec3 intersect = glm::rotateX(glm::rotateY(glm::rotateZ(intersection - center, -rotation.z), -rotation.y), -rotation.x);

        if (glm::length(glm::vec2(intersect.x, intersect.z)) - radius <= FLOAT_ERROR)
            norm = normalize(glm::vec3(intersect.x, 0.f, intersect.z));
        else
            return ZERO_VECTOR;

        norm = glm::rotateZ(glm::rotateY(glm::rotateX(norm, rotation.x), rotation.y), rotation.z);
        return (differenceB) ? -norm : norm;
    }
    void scale(bool enlarge){}
    void move(glm::vec3 move) {}
    void rotate(glm::vec3 rotate) {}
    void breakBoolean(std::vector<Object*> *objectVec, int index) {}
};

// the actual objects
struct Sphere : Object
{
    Sphere()
    {
        radius = BASE_RADIUS;
        center = BASE_CENTER;
        colour = generateRandomVector();
        phong = BASE_PHONG;
    }

    /*
    Finds where the ray enters and exits this sphere, if at all, and adds the points to the rays volume array
    Algorithm from https://www.cs.princeton.edu/courses/archive/fall00/cs426/lectures/raycast/sld013.htm
    */
    void getVolume(Ray *ray)
    {
        Ray tempRay = ray->generateTempRay(center, rotation);
        glm::vec3 relativeCenter = ZERO_VECTOR - tempRay.origin;
        float directionAngle = dot(relativeCenter, tempRay.direction);

        if (directionAngle >= 0)
        {
            float d2 = dot(relativeCenter, relativeCenter) - (directionAngle * directionAngle);
            float r2 = radius * radius;
            if (d2 <= r2)
            {
                float root = sqrt(r2 - d2);
                float scalar1 = directionAngle - root;
                float scalar2 = directionAngle + root;
                ray->pushVolume(glm::min(scalar1, scalar2), glm::max(scalar1, scalar2), ray, this);
            }
        }
    }
    glm::vec3 getNormal(glm::vec3 intersection)
    {
        glm::vec3 intersect = glm::rotateX(glm::rotateY(glm::rotateZ(intersection - center, -rotation.z), -rotation.y), -rotation.x);
        if (abs(glm::length(intersect) - radius) < FLOAT_ERROR)
        {
            glm::vec3 norm = glm::rotateZ(glm::rotateY(glm::rotateX(intersect, rotation.x), rotation.y), rotation.z);
            return (differenceB) ? -glm::normalize(norm) : glm::normalize(norm);
        }
        return ZERO_VECTOR;
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
    void breakBoolean(std::vector<Object*> *objectVec, int index) {}
	void getTree() 
	{
		 
	}
};
struct Cube : Object
{
    // variables
    Plane* planes[6];
    glm::vec3 rotation;
    float radi;

    // constructors
    Cube()
    {
        radius = BASE_RADIUS;
        center = BASE_CENTER;
        colour = generateRandomVector();
        phong = BASE_PHONG;

        planes[0] = new Plane(center + (radius * XAXIS), XAXIS, colour);
        planes[1] = new Plane(center - (radius * XAXIS), -XAXIS, colour);
        planes[2] = new Plane(center + (radius * YAXIS), YAXIS, colour);
        planes[3] = new Plane(center - (radius * YAXIS), -YAXIS, colour);
        planes[4] = new Plane(center + (radius * ZAXIS), ZAXIS, colour);
        planes[5] = new Plane(center - (radius * ZAXIS), -ZAXIS, colour);
    }
    ~Cube()
    {
        for (Plane* plane : planes)
            delete plane;
    };

    void getVolume(Ray *ray)
    {
        Ray tempRay = ray->generateTempRay(center, rotation);

        std::vector<float> points;

        // 6 is the number of planes in planes
        for (Plane* plane : planes)
        {
            float scalar = plane->getIntersection(&tempRay);
            // if parallel, the scalar will be NULL
            if (!scalar) continue;

            // intersect position relative to the cube center
            glm::vec3 intersection = tempRay.applyScalar(scalar);
            glm::vec3 corner = radius * (YAXIS + XAXIS + ZAXIS);

            // if the intersection does not exceed the bounds of the normals
            if (abs(intersection.x) <= abs(corner.x) + FLOAT_ERROR &&
                abs(intersection.y) <= abs(corner.y) + FLOAT_ERROR &&
                abs(intersection.z) <= abs(corner.z) + FLOAT_ERROR)
                points.push_back(scalar);
        }

        if (points.size() < 2) return;

        std::sort(points.begin(), points.end());
        ray->pushVolume(points.front(), points.back(), ray, this);
    }
    glm::vec3 getNormal(glm::vec3 intersection)
    {
        glm::vec3 intersect = glm::rotateX(glm::rotateY(glm::rotateZ(intersection, -rotation.z), -rotation.y), -rotation.x) - center;
        for (Plane *plane : planes)
        {
            if (abs(dot(normalize(intersect - plane->center), plane->normal)) <= FLOAT_ERROR)
            {
                glm::vec3 norm = glm::rotateZ(glm::rotateY(glm::rotateX(plane->normal, rotation.x), rotation.y), rotation.z);
                return (differenceB) ? -norm : norm;
            }
        }
        return ZERO_VECTOR;
    }
    void scale(bool enlarge)
    {
        radius = (enlarge) ? radius + SCALE_CHANGE : glm::max(MIN_SCALE, radius - SCALE_CHANGE);
        for (Plane *plane : planes)
            plane->center = (radius * plane->normal);
    }
    void move(glm::vec3 move) { center += move; }
    void rotate(glm::vec3 rotate) { rotation += rotate; }
    void breakBoolean(std::vector<Object*> *objectVec, int index) {}
};
struct Torus : Object
{
    glm::vec3 rotation;
    float R, r;

    Torus()
    {
        center = BASE_CENTER;
        radius = BASE_RADIUS;
        colour = generateRandomVector();
        phong = BASE_PHONG;
        R = BASE_RADIUS;
        r = BASE_RADIUS * .25f;
    }


    
    //  ALL QUARTIC FUNCTIONS, RAY-TORUS INTERSECTIONS, AND TORUS NORMALS, HAVE BEEN TAKEN FROM THE CODE BASE ON http://www.cosinekitty.com/raytrace/ AND MODIFIED TO FIT THE CODE STYLE USED
    // quartic helper functions
    const double TOLERANCE = 1.0e-8;
    int FilterRealNumbers(int numComplexValues, const std::complex<double> inArray[], double outArray[])
    {
        int numRealValues = 0;
        for (int i = 0; i < numComplexValues; ++i)
            if (fabs(inArray[i].imag()) < TOLERANCE)
                outArray[numRealValues++] = inArray[i].real();
        return numRealValues;
    }
    bool IsZero(std::complex<double> x)
    {
        return (fabs(x.real()) < TOLERANCE) && (fabs(x.imag()) < TOLERANCE);
    }
    std::complex<double> cbrt(std::complex<double> a, int n)
    {
        // this funcion always returns the first complex root
        #define TWOPI (2.0 * 3.141592653589793238462643383279502884)

        double rho = pow(abs(a), 1.0 / 3.0);
        double theta = ((TWOPI * n) + arg(a)) / 3.0;
        return std::complex<double>(rho * cos(theta), rho * sin(theta));
    }
    // Returns n=0..2, the number of distinct real roots found for the equation
    //     ax^2 + bx + c = 0
    // Stores the roots in the first n slots of the array 'roots'.
    int SolveQuadraticEquation(std::complex<double> a, std::complex<double> b, std::complex<double> c, std::complex<double> roots[2])
    {
        if (IsZero(a))
        {
            if (IsZero(b))
                return 0;   // cannot divide by zero, so there is no solution.
            else
            {
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
    // Returns n=0..3, the number of distinct real roots found for the equation
    //     ax^3 + bx^2 + cx + d = 0
    // Stores the roots in the first n slots of the array 'roots'.
    int SolveCubicEquation(std::complex<double> a, std::complex<double> b, std::complex<double> c, std::complex<double> d, std::complex<double> roots[3])
    {
        if (IsZero(a))
            return SolveQuadraticEquation(b, c, d, roots);

        b /= a;
        c /= a;
        d /= a;

        std::complex<double> S = b / 3.0;
        std::complex<double> D = c / 3.0 - S*S;
        std::complex<double> E = S*S*S + (d - S*c) / 2.0;
        std::complex<double> Froot = sqrt(E*E + D*D*D);
        std::complex<double> F = -Froot - E;

        if (IsZero(F))
            F = Froot - E;

        for (int i = 0; i < 3; ++i)
        {
            const std::complex<double> G = cbrt(F, i);
            roots[i] = G - D / G - S;
        }

        return 3;
    }
    // Returns n=0..4, the number of distinct real roots found for the equation
    //     ax^4 + bx^3 + cx^2 + dx + e = 0
    // Stores the roots in the first n slots of the array 'roots'.
    int SolveQuarticEquation(std::complex<double> a, std::complex<double> b, std::complex<double> c, std::complex<double> d, std::complex<double> e, double roots[4])
    {
        std::complex<double> croots[4];

        if (IsZero(a))
        {
            SolveCubicEquation(b, c, d, e, croots);
            return FilterRealNumbers(3, croots, roots);
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

            croots[0] = t + r1;
            croots[1] = t - r1;
            croots[2] = t + r2;
            croots[3] = t - r2;
        }
        else
        {
            std::complex<double> alpha3 = alpha * alpha2;
            std::complex<double> P = -(alpha2 / 12.0 + gamma);
            std::complex<double> Q = -alpha3 / 108.0 + alpha*gamma / 3.0 - beta*beta / 8.0;
            std::complex<double> R = -Q / 2.0 + sqrt(Q*Q / 4.0 + P*P*P / 27.0);
            std::complex<double> U = cbrt(R, 0);
            std::complex<double> y = (-5.0 / 6.0)*alpha + U;

            y -= (IsZero(U)) ? cbrt(Q, 0) : P / (3.0 * U);

            std::complex<double> W = sqrt(alpha + 2.0*y);

            std::complex<double> r1 = sqrt(-(3.0*alpha + 2.0*y + 2.0*beta / W));
            std::complex<double> r2 = sqrt(-(3.0*alpha + 2.0*y - 2.0*beta / W));

            croots[0] = t + (W - r1) / 2.0;
            croots[1] = t + (W + r1) / 2.0;
            croots[2] = t + (-W - r2) / 2.0;
            croots[3] = t + (-W + r2) / 2.0;
        }

        return FilterRealNumbers(4, croots, roots);
    }


    /*
    Finds where the ray interset the trus and adds all volumes to the ray's volume array
    resources for intersection
    http://www.cosinekitty.com/raytrace/chapter13_torus.html

    quartic equation solver has been copied then modified from the code provided on the same website
    */
    void getVolume(Ray *ray)
    {  
        Ray tempRay = ray->generateTempRay(center, rotation);

        // BEGINNING OF MODIFIED SOURCED CODE

        double T = 4.f * R * R;
        double G = T * ((tempRay.direction.x * tempRay.direction.x) + (tempRay.direction.y * tempRay.direction.y));
        double H = 2.f * T * ((tempRay.origin.x * tempRay.direction.x) + (tempRay.origin.y * tempRay.direction.y));
        double I = T * ((tempRay.origin.x * tempRay.origin.x) + (tempRay.origin.y * tempRay.origin.y));
        double J = (tempRay.direction.x * tempRay.direction.x) + (tempRay.direction.y * tempRay.direction.y) + (tempRay.direction.z * tempRay.direction.z);
        double K = 2.f * dot(tempRay.origin, tempRay.direction);
        double L = (tempRay.origin.x * tempRay.origin.x) + (tempRay.origin.y * tempRay.origin.y) + (tempRay.origin.z * tempRay.origin.z) + (R * R) - (r * r);

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
               scalars.push_back(roots[i]);
        }

        // END OF MODIFIED SOURCED CODE


        std::sort(scalars.begin(), scalars.end());
        switch (scalars.size())
        {
        case(2):
            ray->pushVolume(scalars[0], scalars[1], ray, this);
            break;
        case(4):
            ray->pushVolume(scalars[0], scalars[1], ray, this);
            ray->pushVolume(scalars[2], scalars[3], ray, this);
            break;
        default:
            break;
        }
    }
    glm::vec3 getNormal(glm::vec3 intersection)
    {
        glm::vec3 intersect = glm::rotateX(glm::rotateY(glm::rotateZ(intersection - center, -rotation.z), -rotation.y), -rotation.x);


        float implicit = pow(R - sqrt((intersect.x * intersect.x) + (intersect.y * intersect.y)), 2.f) + (intersect.z * intersect.z) - (r * r);

        if (implicit >= FLOAT_ERROR)
            return ZERO_VECTOR;



        float a = 1.f - (R / sqrt((intersect.x * intersect.x) + (intersect.y * intersect.y)));
        glm::vec3 norm = normalize(glm::vec3(a * intersect.x, a * intersect.y, intersect.z));

        norm = glm::rotateZ(glm::rotateY(glm::rotateX(norm, rotation.x), rotation.y), rotation.z);
        return (differenceB) ? -norm : norm;
    }
    void scale(bool enlarge) 
    {
        R += (enlarge) ? SCALE_CHANGE : -SCALE_CHANGE;
        r = R * .25f;
    }
    void move(glm::vec3 move) { center += move; }
    void rotate(glm::vec3 rotate) { rotation += rotate; }
    void breakBoolean(std::vector<Object*> *objectVec, int index) {}
};
struct Cylinder : Object
{
    #define BASE_LENGTH .5f
    float length;
    glm::vec3 rotation = ZERO_VECTOR;
    Plane *topPlane, *bottomPlane;

    Cylinder()
    {
        radius = BASE_RADIUS;
        center = BASE_CENTER;
        colour = generateRandomVector();
        phong = BASE_PHONG;
        length = BASE_LENGTH;

        topPlane = new Plane(center + (length * YAXIS), YAXIS, colour);
        bottomPlane = new Plane(center - (length * YAXIS), -YAXIS, colour);
    }
    ~Cylinder() { delete topPlane; delete bottomPlane; }

    void getVolume(Ray *ray)
    {
        Ray tempRay = ray->generateTempRay(center, rotation);

        std::vector<float> scalars;

        // get and check intersections with planes, if the ray collides with the plane within the bounds of the sphere, add the scalar to the scalar list
        float tempScalar = topPlane->getIntersection(&tempRay);
        if (tempScalar != NULL)
        {
            glm::vec3 intersection = tempRay.applyScalar(tempScalar);
            if (distance(intersection, topPlane->center) - radius <= FLOAT_ERROR)
                scalars.push_back(tempScalar);
        }
        tempScalar = bottomPlane->getIntersection(&tempRay);
        if (tempScalar != NULL)
        {
            glm::vec3 intersection = tempRay.applyScalar(tempScalar);
            if (distance(intersection, bottomPlane->center) - radius <= FLOAT_ERROR)
                scalars.push_back(tempScalar);
        }

        // if the ray collides perfectly with both planes then it is inside the cylinder, and we dont need to check the cylinder
        if (scalars.size() == 2)
        {
            ray->pushVolume(glm::min(scalars[0], scalars[1]), glm::max(scalars[0], scalars[1]), ray, this);
            return;
        }


        // at most 1 plane intersection. Need to check the cylinder now
        float A = (tempRay.direction.x * tempRay.direction.x) + (tempRay.direction.z * tempRay.direction.z);
        float B = 2.f * ((tempRay.origin.x * tempRay.direction.x) + (tempRay.origin.z * tempRay.direction.z));
        float C = (tempRay.origin.x * tempRay.origin.x) + (tempRay.origin.z * tempRay.origin.z) - (radius * radius);

        // now we solve the quadratic equation
        float rootTerm = (B * B) - (4.f * A * C);
        if (rootTerm < 0) return;   // no solutions, misses the cylinder

        rootTerm = sqrt(rootTerm);

        float scalar1 = (-B - rootTerm) / (2.f * A);
        float scalar2 = (-B + rootTerm) / (2.f * A);

        glm::vec3 intersection1 = tempRay.applyScalar(scalar1);
        glm::vec3 intersection2 = tempRay.applyScalar(scalar2);

        if (abs(intersection1.y) - length <= FLOAT_ERROR)
            scalars.push_back(scalar1);
        if (abs(intersection2.y) - length <= FLOAT_ERROR)
            scalars.push_back(scalar2);

        if (scalars.size() < 2) return;

        std::sort(scalars.begin(), scalars.end());
        ray->pushVolume(scalars.front(), scalars.back(), ray, this);
    }
    glm::vec3 getNormal(glm::vec3 intersection)
    {
        glm::vec3 norm(0.f);
        glm::vec3 intersect= glm::rotateX(glm::rotateY(glm::rotateZ(intersection - center, -rotation.z), -rotation.y), -rotation.x);


        // top plane
        if (abs(intersect.y - topPlane->center.y) <= FLOAT_ERROR)
            norm = topPlane->normal;
        // bottom plane
        else if (abs(intersect.y - bottomPlane->center.y) <= FLOAT_ERROR)
            norm = bottomPlane->normal;
        //cylinder
        else if (glm::length(glm::vec2(intersect.x, intersect.z)) - radius <= FLOAT_ERROR)
            norm = normalize(glm::vec3(intersect.x, 0.f, intersect.z));
        else
            return ZERO_VECTOR;


        norm = glm::rotateZ(glm::rotateY(glm::rotateX(norm, rotation.x), rotation.y), rotation.z);
        return (differenceB) ? -norm : norm;
    }
    void scale(bool enlarge)
    {
        if (enlarge)
        {
            radius += SCALE_CHANGE;
            length += SCALE_CHANGE;
            topPlane->center.y += SCALE_CHANGE;
            bottomPlane->center.y += -SCALE_CHANGE;
        } 
        else
        {
            radius = glm::max(MIN_SCALE, radius - SCALE_CHANGE);
            length = glm::max(MIN_SCALE, length - SCALE_CHANGE);
            topPlane->center.y = ((topPlane->center.y <= .1f) ? topPlane->center.y : topPlane->center.y - SCALE_CHANGE);
            bottomPlane->center.y = ((bottomPlane->center.y >= -.1f) ? bottomPlane->center.y : bottomPlane->center.y + SCALE_CHANGE);
        }
    }
    void move(glm::vec3 move) { center += move; }
    void rotate(glm::vec3 rotate) { rotation += rotate; }
    void breakBoolean(std::vector<Object*> *objectVec, int index) {}
};

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// the combinations

struct Union : Object
{
    Object *leftChild, *rightChild;
    Union(Object *obj1, Object *obj2) : leftChild(obj1), rightChild(obj2)
    {
        phong = (obj1->phong + obj2->phong) / 2.f;
        colour = (obj1->colour + obj2->colour) / 2.f;
        center = (obj1->center + obj2->center) / 2.f;
        selected = false;
    }

    void getVolume(Ray *ray)
    {
        Ray tempRay1(ray->origin, ray->direction);
        Ray tempRay2(ray->origin, ray->direction);

        // get all the volumes
        leftChild->getVolume(&tempRay1);
        rightChild->getVolume(&tempRay2);

        // create seperate lists of volumes for easier access
        std::vector<Volume> v1 = tempRay1.getVolume();
        std::vector<Volume> v2 = tempRay2.getVolume();

        // ray hits neither object
        if (v1.size() == 0 && v2.size() == 0) return;

        // handle the cases where one of the lists is empty
        if (v1.size() == 0)
        {
            for (Volume vol : v2)
                ray->pushVolumeBoolean(vol.entrance, vol.exit, vol.entranceNormal, vol.exitNormal, this);
            return;
        }
        if (v2.size() == 0)
        {
            for (Volume vol : v1)
                ray->pushVolumeBoolean(vol.entrance, vol.exit, vol.entranceNormal, vol.exitNormal, this);
            return;
        }


        // index trackers for each list
        int v1Index = 0;
        int v2Index = 0;

        // temporary storge for enterences and exits
        // initialize the floats to smallest possible
        float tempEntr = -FLT_MAX;
        float tempExit = -FLT_MAX;
        Object *tempObjEntr = NULL;
        Object *tempObjExit = NULL;

        // while both lists still have unvisited elements
        while (true)
        {
            if (v1Index >= v1.size() && v2Index >= v2.size())
                break;

            // variables for easier access of volumes, initialize to zero as there is no base construct
            Volume vol1 = v1[0], vol2 = v2[0];

            // base loop condition. both lists still contatain volumes
            if (v1Index < v1.size() && v2Index < v2.size())
            {
                vol1 = v1[v1Index];
                vol2 = v2[v2Index];
            }
            // if v1 is spent
            else if (v1Index >= v1.size())
            {
                vol1 = Volume(FLT_MAX, FLT_MAX, ZERO_VECTOR, ZERO_VECTOR, this);
                vol2 = v2[v2Index];
            }
            // if v2 is spent
            else if (v2Index >= v2.size())
            {
                vol1 = v1[v1Index];
                vol2 = Volume(FLT_MAX, FLT_MAX, ZERO_VECTOR, ZERO_VECTOR, this);
            }






            // if vol1 is closer than vol2
            if (vol1.entrance < vol2.entrance)
            {
                // vol1 and vol2 are disjoint
                if (vol1.exit < vol2.entrance)
                {
                    // old volume ends before vol1 begins
                    if (tempExit < vol1.entrance)
                    {
                        if (tempExit != -FLT_MAX)
                            ray->pushVolumeBoolean(tempEntr, tempExit, tempObjEntr->getNormal(ray->applyScalar(tempEntr)), tempObjExit->getNormal(ray->applyScalar(tempExit)), this);

                        tempEntr = vol1.entrance;
                        tempExit = vol1.exit;
                        tempObjEntr = vol1.object;
                        tempObjExit = vol1.object;

                        v1Index++;
                    }
                    // old volume ends before vol1 ends
                    else if (tempExit < vol1.exit)
                    {
                        tempExit = vol1.exit;
                        tempObjExit = vol1.object;

                        v1Index++;
                    }
                    // old volume ends before vol2 starts. vol1 is spent
                    else if (tempExit < vol2.entrance)
                    {
                        v1Index++;
                    }
                    // old volume ends before vol2 ends. vol1 is spent. extend to vol2. vol2 then spent
                    else if (tempExit < vol2.exit)
                    {
                        tempExit = vol2.exit;
                        tempObjExit = vol2.object;

                        v1Index++;
                        v2Index++;
                    }
                    // old volume extends past vol2's exit. both vol1 and vol2 are spent
                    else
                    {
                        v1Index++;
                        v2Index++;
                    }
                }
                // vol1 and vol2 partially overlap
                else if (vol1.exit < vol2.exit)
                {
                    // old volume ends before vol1 starts
                    if (tempExit < vol1.entrance)
                    {
                        if (tempExit != -FLT_MAX)
                            ray->pushVolumeBoolean(tempEntr, tempExit, tempObjEntr->getNormal(ray->applyScalar(tempEntr)), tempObjExit->getNormal(ray->applyScalar(tempExit)), this);

                        tempEntr = vol1.entrance;
                        tempExit = vol2.exit;
                        tempObjEntr = vol1.object;
                        tempObjExit = vol2.object;

                        v1Index++;
                        v2Index++;
                    }
                    // old volume ends before vol2 begins
                    else if (tempExit < vol2.entrance)
                    {
                        tempExit = vol2.exit;
                        tempObjExit = vol2.object;

                        v1Index++;
                        v2Index++;
                    }
                    // old volume ends before vol1 ends
                    else if (tempExit < vol1.exit)
                    {
                        tempExit = vol2.exit;
                        tempObjExit = vol2.object;

                        v1Index++;
                        v2Index++;
                    }
                    // old volume ends before vol2 ends
                    else if (tempExit < vol2.exit)
                    {
                        tempExit = vol2.exit;
                        tempObjExit = vol2.object;

                        v1Index++;
                        v2Index++;
                    }
                    // old volume extends past vol2
                    else
                    {
                        v1Index++;
                        v2Index++;
                    }
                }
                // vol2 encompasses vol2
                else
                {
                    // old volume ends before vol1 begins
                    if (tempExit < vol1.entrance)
                    {
                        if (tempExit != -FLT_MAX)
                            ray->pushVolumeBoolean(tempEntr, tempExit, tempObjEntr->getNormal(ray->applyScalar(tempEntr)), tempObjExit->getNormal(ray->applyScalar(tempExit)), this);

                        tempEntr = vol1.entrance;
                        tempExit = vol1.exit;
                        tempObjEntr = vol1.object;
                        tempObjExit = vol1.object;

                        v1Index++;
                        v2Index++;
                    }
                    // old volume ends before vol2 begins
                    else if (tempExit < vol2.entrance)
                    {
                        tempExit = vol1.exit;
                        tempObjExit = vol1.object;

                        v1Index++;
                        v2Index++;
                    }
                    // old volume ends before vol2 ends
                    else if (tempExit < vol2.exit)
                    {
                        tempExit = vol1.exit;
                        tempObjExit = vol1.object;

                        v1Index++;
                        v2Index++;
                    }
                    // old volume ends before vol1 ends
                    else if (tempExit < vol1.exit)
                    {
                        tempExit = vol1.exit;
                        tempObjExit = vol1.object;

                        v1Index++;
                        v2Index++;
                    }
                    // old volume extends past vol1's end
                    else
                    {
                        v1Index++;
                        v2Index++;
                    }
                }
            }
            // if vol2 is closer than vol1
            else
            {
                // vol1 and vol2 are disjoint
                if (vol2.exit < vol1.entrance)
                {
                    // old volume ends before vol2 begins
                    if (tempExit < vol2.entrance)
                    {
                        if (tempExit != -FLT_MAX)
                            ray->pushVolumeBoolean(tempEntr, tempExit, tempObjEntr->getNormal(ray->applyScalar(tempEntr)), tempObjExit->getNormal(ray->applyScalar(tempExit)), this);

                        tempEntr = vol2.entrance;
                        tempExit = vol2.exit;
                        tempObjEntr = vol2.object;
                        tempObjExit = vol2.object;

                        v2Index++;
                    }
                    // old volume ends before vol2 ends
                    else if (tempExit < vol2.exit)
                    {
                        tempExit = vol2.exit;
                        tempObjExit = vol2.object;

                        v2Index++;
                    }
                    // old volume ends before vol1 begins
                    else if (tempExit < vol1.entrance)
                    {
                        v2Index++;
                    }
                    // old volume ends before vol1 ends
                    else if (tempExit < vol1.exit)
                    {
                        tempExit = vol1.exit;
                        tempObjExit = vol1.object;

                        v1Index++;
                        v2Index++;
                    }
                    // old volume extends past vol1
                    else
                    {
                        v1Index++;
                        v2Index++;
                    }
                }
                // vol1 and vol2 partially overlap
                else if (vol2.exit < vol1.exit)
                {
                    // old volume ends before vol2 begins
                    if (tempExit < vol2.entrance)
                    {
                        if (tempExit != -FLT_MAX)
                            ray->pushVolumeBoolean(tempEntr, tempExit, tempObjEntr->getNormal(ray->applyScalar(tempEntr)), tempObjExit->getNormal(ray->applyScalar(tempExit)), this);

                        tempEntr = vol2.entrance;
                        tempExit = vol1.exit;
                        tempObjEntr = vol2.object;
                        tempObjExit = vol1.object;

                        v1Index++;
                        v2Index++;
                    }
                    // old volume ends before vol1 begins
                    else if (tempExit < vol1.entrance)
                    {
                        tempExit = vol1.exit;
                        tempObjExit = vol1.object;

                        v1Index++;
                        v2Index++;
                    }
                    // old volume ends before vol2 ends
                    else if (tempExit < vol2.exit)
                    {
                        tempExit = vol1.exit;
                        tempObjExit = vol1.object;

                        v1Index++;
                        v2Index++;
                    }
                    // old volume ends before vol1 ends
                    else if (tempExit < vol1.exit)
                    {
                        tempExit = vol1.exit;
                        tempObjExit = vol1.object;

                        v1Index++;
                        v2Index++;
                    }
                    // old volume extends past vol1
                    else
                    {
                        v1Index++;
                        v2Index++;
                    }
                }
                // vol2 encompasses vol1
                else
                {
                    // old volume ends before vol2 begins
                    if (tempExit < vol2.entrance)
                    {
                        if (tempExit != -FLT_MAX)
                            ray->pushVolumeBoolean(tempEntr, tempExit, tempObjEntr->getNormal(ray->applyScalar(tempEntr)), tempObjExit->getNormal(ray->applyScalar(tempExit)), this);

                        tempEntr = vol2.entrance;
                        tempExit = vol2.exit;
                        tempObjEntr = vol2.object;
                        tempObjExit = vol2.object;

                        v1Index++;
                        v2Index++;
                    }
                    // old volume ends before vol1 begins
                    else if (tempExit < vol1.entrance)
                    {
                        tempExit = vol2.exit;
                        tempObjExit = vol2.object;

                        v1Index++;
                        v2Index++;
                    }
                    // old volume ends before vol1 ends
                    else if (tempExit < vol1.exit)
                    {
                        tempExit = vol2.exit;
                        tempObjExit = vol2.object;

                        v1Index++;
                        v2Index++;
                    }
                    // old volume ends before vol2 ends
                    else if (tempExit < vol2.exit)
                    {
                        tempExit = vol2.exit;
                        tempObjExit = vol2.object;

                        v1Index++;
                        v2Index++;
                    }
                    // old volume extends past vol2
                    else
                    {
                        v1Index++;
                        v2Index++;
                    }
                }
            }

        }

        // push the remaining volume
        if (tempExit != -FLT_MAX)
            ray->pushVolumeBoolean(tempEntr, tempExit, tempObjEntr->getNormal(ray->applyScalar(tempEntr)), tempObjExit->getNormal(ray->applyScalar(tempExit)), this);
    }
    glm::vec3 getNormal(glm::vec3 intersection) 
    {
        glm::vec3 lNormal = leftChild->getNormal(intersection);
        if (lNormal != ZERO_VECTOR)
            return (differenceB) ? -lNormal : lNormal;

        glm::vec3 rNormal = rightChild->getNormal(intersection);
        if (rNormal != ZERO_VECTOR)
            return (differenceB) ? -rNormal : rNormal;

        return ZERO_VECTOR;
    }
    void scale(bool enlarge){}
    void move(glm::vec3 move)
    {
        leftChild->move(move);
        rightChild->move(move);
    }
    void rotate(glm::vec3 rotate){}
    void select() { selected = true; }
    void deselect() { selected = false; }
    void breakBoolean(std::vector<Object*> *objectVec, int index)
    {
        objectVec->erase(objectVec->begin() + index);
        objectVec->push_back(leftChild);
        objectVec->push_back(rightChild);
        delete this;
    }
};
struct Intersection : Object
{
    Object *leftChild, *rightChild;
    Intersection(Object *obj1, Object *obj2) : leftChild(obj1), rightChild(obj2)
    {
        phong = (obj1->phong + obj2->phong) / 2.f;
        colour = (obj1->colour + obj2->colour) / 2.f;
        center = (obj1->center + obj2->center) / 2.f;
        selected = false;
    }

    void getVolume(Ray *ray)
    {
        Ray tempRay1(ray->origin, ray->direction);
        Ray tempRay2(ray->origin, ray->direction);

        // get all the volumes
        leftChild->getVolume(&tempRay1);
        rightChild->getVolume(&tempRay2);

        // create seperate lists of volumes for easier access
        std::vector<Volume> v1 = tempRay1.getVolume();
        std::vector<Volume> v2 = tempRay2.getVolume();

        // ray hits neither, or only 1 object then intersection does not care
        if (v1.size() == 0 || v2.size() == 0) return;

        // intex trackers for each list
        int v1Index = 0;
        int v2Index = 0;

        // temporary storge for enterences and exits, initialize the floats to smallest possible
        float tempEntr = -FLT_MAX;
        float tempExit = -FLT_MAX;
        Object *tempObjEntr = NULL;
        Object *tempObjExit = NULL;

        while (v1Index != v1.size() && v2Index != v2.size())
        {
            // define as new volume for easier access
            Volume vol1 = v1[v1Index], vol2 = v2[v2Index];
            
            // vol1 is closer than vol2
            if (vol1.entrance < vol2.entrance)
            {
                // vol 1 and vol2 are disjoint. the disjointness allows us to push the old volume, reset the defaults, and increment v1
                if (vol1.exit < vol2.entrance)
                {
                    // old volume ends before vol2 begins. before a possible intersection. vol1 is spent. vol2 has potential
                    if (tempExit < vol2.entrance)
                    {
                        if (tempExit != -FLT_MAX)
                            ray->pushVolumeBoolean(tempEntr, tempExit, tempObjEntr->getNormal(ray->applyScalar(tempEntr)), tempObjExit->getNormal(ray->applyScalar(tempExit)), this);

                        tempEntr = -FLT_MAX;
                        tempExit = -FLT_MAX;
                        tempObjEntr = NULL;
                        tempObjExit = NULL;

                        v1Index++;
                    }
                    // old volume ends before vol2 ends. vol1 is spent. vol2 has potential
                    else if (tempExit < vol2.exit)
                    {
                        v1Index++;
                    }
                    // old volume extends post vol2. both vol1 and vol2 are spent.
                    else
                    {
                        v1Index++;
                        v2Index++;
                    }
                }
                // vol1 and vol2 partially overlap
                else if (vol1.exit < vol2.exit)
                {
                    // old volume ends before vol1 starts
                    if (tempExit < vol1.entrance)
                    {
                        if (tempExit != -FLT_MAX)
                            ray->pushVolumeBoolean(tempEntr, tempExit, tempObjEntr->getNormal(ray->applyScalar(tempEntr)), tempObjExit->getNormal(ray->applyScalar(tempExit)), this);

                        tempEntr = vol2.entrance;
                        tempExit = vol1.exit;
                        tempObjEntr = vol2.object;
                        tempObjExit = vol1.object;

                        v1Index++;
                    }
                    // old volume ends before vol2 starts
                    else if (tempExit < vol2.entrance)
                    {
                        ray->pushVolumeBoolean(tempEntr, tempExit, tempObjEntr->getNormal(ray->applyScalar(tempEntr)), tempObjExit->getNormal(ray->applyScalar(tempExit)), this);

                        tempEntr = vol2.entrance;
                        tempExit = vol1.exit;
                        tempObjEntr = vol2.object;
                        tempObjExit = vol1.object;

                        v1Index++;
                    }
                    // old volume ends before vol1 ends
                    else if (tempExit < vol1.exit)
                    {
                        ray->pushVolumeBoolean(tempEntr, vol1.exit, tempObjEntr->getNormal(ray->applyScalar(tempEntr)), vol1.exitNormal, this);

                        tempEntr = -FLT_MAX;
                        tempExit = -FLT_MAX;
                        tempObjEntr = NULL;
                        tempObjExit = NULL;

                        v1Index++;
                    }
                    // old volume ends before vol2 ends
                    else if (tempExit < vol2.exit)
                    {
                        // no intersection is added, but there is room for growth. vol1 passed
                        v1Index++;
                    }
                    // old volume extends past vol2's end
                    else 
                    {
                        // no intersections added. no possibility for vol1 or vol2 to add an intersection
                        v1Index++;
                        v2Index++;
                    }
                }
                // vol1 completely overlaps vol2
                else
                {
                    // old volume ends before vol1 begins. v2 gets encompased vy vol1. start new temp with the intersection
                    if (tempExit < vol1.entrance)
                    {
                        if (tempExit != -FLT_MAX)
                            ray->pushVolumeBoolean(tempEntr, tempExit, tempObjEntr->getNormal(ray->applyScalar(tempEntr)), tempObjExit->getNormal(ray->applyScalar(tempExit)), this);

                        tempEntr = vol2.entrance;
                        tempExit = vol2.exit;
                        tempObjEntr = vol2.object;
                        tempObjExit = vol2.object;

                        v2Index++;
                    }
                    // old volume ends before vol2 begins.
                    else if (tempExit < vol2.entrance)
                    {
                        ray->pushVolumeBoolean(tempEntr, tempExit, tempObjEntr->getNormal(ray->applyScalar(tempEntr)), tempObjExit->getNormal(ray->applyScalar(tempExit)), this);

                        tempEntr = vol2.entrance;
                        tempExit = vol2.exit;
                        tempObjEntr = vol2.object;
                        tempObjExit = vol2.object;

                        v2Index++;
                    }
                    // old volume ends before vol2 ends. extend the old volume, vol2 is spent
                    else if (tempExit < vol2.exit)
                    {
                        tempExit = vol2.exit;
                        tempObjExit = vol2.object;

                        v2Index++;
                    }
                    // old volume ends before vol1 ends. vol1 U vol2 add nothing. vol2 is spent. vol1 still has potential
                    else if (tempExit < vol1.exit)
                    {
                        v2Index++;
                    }
                    // old volume ends after vol1 ends. nothing added. both spent
                    else
                    {
                        v1Index++;
                        v2Index++;
                    }
                }
            }
            // vol2 is closer than vol1
            else
            {
                // vol2 and vol1 are disjoint.
                if (vol2.exit < vol1.entrance)
                {
                    // no intersections can be added.old volume is done.vol2 is done
                    if (tempExit < vol1.entrance)
                    {
                        if (tempExit != -FLT_MAX)
                            ray->pushVolumeBoolean(tempEntr, tempExit, tempObjEntr->getNormal(ray->applyScalar(tempEntr)), tempObjExit->getNormal(ray->applyScalar(tempExit)), this);

                        tempEntr = -FLT_MAX;
                        tempExit = -FLT_MAX;
                        tempObjEntr = NULL;
                        tempObjExit = NULL;

                        v2Index++;
                    }
                    // old volume extends past vol1's start so vol1 still has potential. old volume can still be added to. vol2 spent
                    else if (tempExit < vol1.exit)
                    {
                        v2Index++;
                    }
                    // old volume extends past vol1. both vol1 and 2 are spent
                    else
                    {
                        v1Index++;
                        v2Index++;
                    }
                }
                // vol2 and vol1 partially overlap
                else if (vol2.exit < vol1.exit)
                {
                    // old volume ends before vol2 begins
                    if (tempExit < vol2.entrance)
                    {
                        if (tempExit != -FLT_MAX)
                            ray->pushVolumeBoolean(tempEntr, tempExit, tempObjEntr->getNormal(ray->applyScalar(tempEntr)), tempObjExit->getNormal(ray->applyScalar(tempExit)), this);

                        tempEntr = vol1.entrance;
                        tempExit = vol2.exit;
                        tempObjEntr = vol1.object;
                        tempObjExit = vol2.object;

                        v2Index++;
                    }
                    // old volume ends bofore vol1 begins
                    else if (tempExit < vol1.entrance)
                    {
                        if (tempExit != -FLT_MAX)
                            ray->pushVolumeBoolean(tempEntr, tempExit, tempObjEntr->getNormal(ray->applyScalar(tempEntr)), tempObjExit->getNormal(ray->applyScalar(tempExit)), this);

                        tempEntr = vol1.entrance;
                        tempExit = vol2.exit;
                        tempObjEntr = vol1.object;
                        tempObjExit = vol2.object;

                        v2Index++;
                    }
                    // old volume ends before vol2 ends
                    else if (tempExit < vol2.exit)
                    {
                        tempExit = vol2.exit;
                        tempObjExit = vol2.object;

                        v2Index++;
                    }
                    // old volume ends before vol1 ends. potential intersection is encompassed. vol2 is spent. vol1 has potential
                    else if (tempExit < vol1.exit)
                    {
                        v2Index++;
                    }
                    // old volume ends after vol1 ends. any intersections are encompassed. both vol1 and vol2 are spent
                    else
                    {
                        v1Index++;
                        v2Index++;
                    }
                }
                // vol2 encompasses vol1
                else
                {
                    // old volume ends before vol2 starts. vol1 is spent. vol2 has potential
                    if (tempExit < vol2.entrance)
                    {
                        if (tempExit != -FLT_MAX)
                            ray->pushVolumeBoolean(tempEntr, tempExit, tempObjEntr->getNormal(ray->applyScalar(tempEntr)), tempObjExit->getNormal(ray->applyScalar(tempExit)), this);

                        tempEntr = vol1.entrance;
                        tempExit = vol1.exit;
                        tempObjEntr = vol1.object;
                        tempObjExit = vol1.object;

                        v1Index++;
                    }
                    // old volume ends before vol1 begins
                    else if (tempExit < vol1.entrance)
                    {
                        ray->pushVolumeBoolean(tempEntr, tempExit, tempObjEntr->getNormal(ray->applyScalar(tempEntr)), tempObjExit->getNormal(ray->applyScalar(tempExit)), this);

                        tempEntr = vol1.entrance;
                        tempExit = vol1.exit;
                        tempObjEntr = vol1.object;
                        tempObjExit = vol1.object;

                        v1Index++;
                    }
                    // old volume ends before vol1 ends
                    else if (tempExit < vol1.exit)
                    {
                        tempExit = vol1.exit;
                        tempObjExit = vol1.object;

                        v1Index++;
                    }
                    // old volume ends before vol2 ends. intersection is encompassed. vol1 is spent. vol2 has potential
                    else if (tempExit < vol2.exit)
                    {
                        v1Index++;
                    }
                    // old volume extends past vol2's exit. intersection encompassed. both vols spent
                    else
                    {
                        v1Index++;
                        v2Index++;
                    }
                }
            }
        }

        // push the remaining volume
        if (tempExit != -FLT_MAX)
            ray->pushVolumeBoolean(tempEntr, tempExit, tempObjEntr->getNormal(ray->applyScalar(tempEntr)), tempObjExit->getNormal(ray->applyScalar(tempExit)), this);
    }
    glm::vec3 getNormal(glm::vec3 intersection)
    {
        glm::vec3 lNormal = leftChild->getNormal(intersection);
        if (lNormal != ZERO_VECTOR)
            return (differenceB) ? -lNormal : lNormal;

        glm::vec3 rNormal = rightChild->getNormal(intersection);
        if (rNormal != ZERO_VECTOR)
            return (differenceB) ? -rNormal : rNormal;

        return ZERO_VECTOR;
    }
    void scale(bool enlarge){}
    void move(glm::vec3 move)
    {
        leftChild->move(move);
        rightChild->move(move);
    }
    void rotate(glm::vec3 rotate){}
    void breakBoolean(std::vector<Object*> *objectVec, int index)
    {
        objectVec->erase(objectVec->begin() + index);
        objectVec->push_back(leftChild);
        objectVec->push_back(rightChild);
        delete this;
    }
};
struct Difference : Object
{
    Object *leftChild, *rightChild;
    Difference(Object *obj1, Object *obj2) : leftChild(obj1), rightChild(obj2)
    {
        phong = (obj1->phong + obj2->phong) / 2.f;
        colour = (obj1->colour + obj2->colour) / 2.f;
        center = (obj1->center + obj2->center) / 2.f;
        selected = false;

        rightChild->differenceB = true;
    }

    //%%%%%%%%%%%%%%%%%%%%%
    // difference = A-B where leftChild=A and rightChild=B
    void getVolume(Ray *ray)
    {
        Ray tempRay1(ray->origin, ray->direction);
        Ray tempRay2(ray->origin, ray->direction);

        // get all the volumes
        leftChild->getVolume(&tempRay1);
        rightChild->getVolume(&tempRay2);

        // create seperate lists of volumes for easier access
        std::vector<Volume> v1 = tempRay1.getVolume();
        std::vector<Volume> v2 = tempRay2.getVolume();

        // ray doesnt hit A. there is no volume
        if (v1.size() == 0) return;

        // B is empty. just keep A's volume
        if (v2.size() == 0)
        {
            for (Volume vol : v1)
                ray->pushVolumeBoolean(vol.entrance, vol.exit, vol.entranceNormal, vol.exitNormal, this);
            return;
        }


        // index trackers for each list
        int v1Index = 0;
        int v2Index = 0;

        // temporary storge for enterences and exits
        // initialize the floats to smallest possible
        float tempEntr = -FLT_MAX;
        float tempExit = -FLT_MAX;
        Object *tempObjEntr = NULL;
        Object *tempObjExit = NULL;

        // vol1 = A,vol2 = B. A-B
        while (true)
        {
            // defaults as volume has no default constructor
            Volume vol1 = v1[0], vol2 = v2[0];

            if (v1Index >= v1.size())
                break;
            else if (v2Index >= v2.size())
            {
                vol1 = v1[v1Index];
                vol2 = Volume(FLT_MAX, FLT_MAX, ZERO_VECTOR, ZERO_VECTOR, this);
            }
            else
            {
                vol1 = v1[v1Index];
                vol2 = v2[v2Index];
            }

            // vol1 is closer than vol2
            if (vol1.entrance < vol2.entrance)
            {
                // vol1 and vol2 are disjoint
                if (vol1.exit < vol2.entrance)
                {
                    // old volume ends before vol1 begins
                    if (tempExit < vol1.entrance)
                    {
                        if (tempExit != -FLT_MAX)
                            ray->pushVolumeBoolean(tempEntr, tempExit, tempObjEntr->getNormal(ray->applyScalar(tempEntr)), tempObjExit->getNormal(ray->applyScalar(tempExit)), this);

                        tempEntr = vol1.entrance;
                        tempExit = vol1.exit;
                        tempObjEntr = vol1.object;
                        tempObjExit = vol1.object;

                        v1Index++;
                    }
                    // old volume ends before vol1 ends
                    else if (tempExit < vol1.exit)
                    {
                        tempExit = vol1.exit;
                        tempObjExit = vol1.object;

                        v1Index++;
                    }
                    // old volume ends before vol2 starts
                    else if (tempExit < vol2.entrance)
                    {
                        v1Index++;
                    }
                    // old volume ends before vol2 ends
                    else if (tempExit < vol2.exit)
                    {
                        if (tempExit != -FLT_MAX)
                            ray->pushVolumeBoolean(tempEntr, vol2.entrance, tempObjEntr->getNormal(ray->applyScalar(tempEntr)), vol2.entranceNormal, this);

                        tempEntr = -FLT_MAX;
                        tempExit = -FLT_MAX;
                        tempObjEntr = NULL;
                        tempObjExit = NULL;

                        v1Index++;
                    }
                    // old volume extends past vol2
                    else
                    {
                        if (tempExit != -FLT_MAX)
                            ray->pushVolumeBoolean(tempEntr, vol2.entrance, tempObjEntr->getNormal(ray->applyScalar(tempEntr)), vol2.entranceNormal, this);

                        tempEntr = vol2.exit;
                        tempObjEntr = vol2.object;
                        v1Index++;
                    }
                }
                // vol1 and vol2 partially overlap
                else if (vol1.exit < vol2.exit)
                {
                    // old volume ends before vol1 begins
                    if (tempExit < vol1.entrance)
                    {
                        if (tempExit != -FLT_MAX)
                            ray->pushVolumeBoolean(tempEntr, tempExit, tempObjEntr->getNormal(ray->applyScalar(tempEntr)), tempObjExit->getNormal(ray->applyScalar(tempExit)), this);

                        tempEntr = vol1.entrance;
                        tempExit = vol2.entrance;
                        tempObjEntr = vol1.object;
                        tempObjExit = vol2.object;

                        v1Index++;
                    }
                    // old volume ends before vol2 begins
                    else if (tempExit < vol2.entrance)
                    {
                        if (tempExit != -FLT_MAX)
                            ray->pushVolumeBoolean(tempEntr, vol2.entrance, tempObjEntr->getNormal(ray->applyScalar(tempEntr)), vol2.entranceNormal, this);

                        tempEntr = -FLT_MAX;
                        tempExit = -FLT_MAX;
                        tempObjEntr = NULL;
                        tempObjExit = NULL;

                        v1Index++;
                    }
                    // old volume ends before vol1 ends
                    else if (tempExit < vol1.exit)
                    {
                        if (tempExit != -FLT_MAX)
                            ray->pushVolumeBoolean(tempEntr, vol2.entrance, tempObjEntr->getNormal(ray->applyScalar(tempEntr)), vol2.entranceNormal, this);

                        tempEntr = -FLT_MAX;
                        tempExit = -FLT_MAX;
                        tempObjEntr = NULL;
                        tempObjExit = NULL;

                        v1Index++;
                    }
                    // old volume ends  before vol2 ends
                    else if (tempExit < vol2.exit)
                    {
                        if (tempExit != -FLT_MAX)
                            ray->pushVolumeBoolean(tempEntr, vol2.entrance, tempObjEntr->getNormal(ray->applyScalar(tempEntr)), vol2.entranceNormal, this);

                        tempEntr = -FLT_MAX;
                        tempExit = -FLT_MAX;
                        tempObjEntr = NULL;
                        tempObjExit = NULL;

                        v1Index++;
                    }
                    // old volume extends past vol2
                    else
                    {
                        if (tempExit != -FLT_MAX)
                            ray->pushVolumeBoolean(tempEntr, vol2.entrance, tempObjEntr->getNormal(ray->applyScalar(tempEntr)), vol2.entranceNormal, this);

                        tempEntr = vol2.exit;
                        tempObjEntr = vol2.object;
                        v1Index++;
                    }
                }
                // vol1 encompasses vol2
                else
                {
                    // old volume ends before vol1 begins
                    if (tempExit < vol1.entrance)
                    {
                        if (tempExit != -FLT_MAX)
                            ray->pushVolumeBoolean(tempEntr, tempExit, tempObjEntr->getNormal(ray->applyScalar(tempEntr)), tempObjExit->getNormal(ray->applyScalar(tempExit)), this);

                        ray->pushVolumeBoolean(vol1.entrance, vol2.entrance, vol1.entranceNormal, vol2.entranceNormal, this);

                        tempEntr = vol2.exit;
                        tempExit = vol1.exit;
                        tempObjEntr = vol2.object;
                        tempObjExit = vol1.object;

                        v1Index++;
                        v2Index++;
                    }
                    // old volume ends before vol2 begins
                    else if (tempExit < vol2.entrance)
                    {
                        if (tempExit != -FLT_MAX)
                            ray->pushVolumeBoolean(tempEntr, vol2.entrance, tempObjEntr->getNormal(ray->applyScalar(tempEntr)), vol2.entranceNormal, this);

                        tempEntr = vol2.exit;
                        tempExit = vol1.exit;
                        tempObjEntr = vol2.object;
                        tempObjExit = vol1.object;

                        v1Index++;
                        v2Index++;
                    }
                    // old volume ends before vol2 ends
                    else if (tempExit < vol2.exit)
                    {
                        if (tempExit != -FLT_MAX)
                            ray->pushVolumeBoolean(tempEntr, vol2.entrance, tempObjEntr->getNormal(ray->applyScalar(tempEntr)), vol2.entranceNormal, this);

                        tempEntr = vol2.exit;
                        tempExit = vol1.exit;
                        tempObjEntr = vol2.object;
                        tempObjExit = vol1.object;

                        v1Index++;
                        v2Index++;
                    }
                    // old volume ends before vol1 ends
                    else if (tempExit < vol1.exit)
                    {
                        if (tempExit != -FLT_MAX)
                            ray->pushVolumeBoolean(tempEntr, vol2.entrance, tempObjEntr->getNormal(ray->applyScalar(tempEntr)), vol2.entranceNormal, this);

                        tempEntr = vol2.exit;
                        tempExit = vol1.exit;
                        tempObjEntr = vol2.object;
                        tempObjExit = vol1.object;

                        v1Index++;
                        v2Index++;
                    }
                    // old volume extends past vol1
                    else
                    {
                        if (tempExit != -FLT_MAX)
                            ray->pushVolumeBoolean(tempEntr, vol2.entrance, tempObjEntr->getNormal(ray->applyScalar(tempEntr)), vol2.entranceNormal, this);

                        tempEntr = vol2.exit;
                        tempObjEntr = vol2.object;

                        v1Index++;
                        v2Index++;
                    }
                }
            }
            // vol2 is closer than vol1
            else
            {
                // vol1 and vol2 are disjoint
                if (vol2.exit < vol1.entrance)
                {
                    // old volume ends before vol2 begins
                    if (tempExit < vol2.entrance)   //CHECK
                    {
                        if (tempExit != -FLT_MAX)
                            ray->pushVolumeBoolean(tempEntr, tempExit, tempObjEntr->getNormal(ray->applyScalar(tempEntr)), tempObjExit->getNormal(ray->applyScalar(tempExit)), this);

                        tempEntr = -FLT_MAX;
                        tempExit = -FLT_MAX;
                        tempObjEntr = NULL;
                        tempObjExit = NULL;

                        v2Index++;
                    }
                    // old volume ends before vol2 ends
                    else if (tempExit < vol2.exit)
                    {
                        if (tempExit != -FLT_MAX)
                            ray->pushVolumeBoolean(tempEntr, vol2.entrance, tempObjEntr->getNormal(ray->applyScalar(tempEntr)), vol2.entranceNormal, this);

                        tempEntr = -FLT_MAX;
                        tempExit = -FLT_MAX;
                        tempObjEntr = NULL;
                        tempObjExit = NULL;

                        v2Index++;
                    }
                    // old volume ends before vol1 begins
                    else if (tempExit < vol1.entrance)
                    {
                        if (tempExit != -FLT_MAX)
                            ray->pushVolumeBoolean(tempEntr, vol2.entrance, tempObjEntr->getNormal(ray->applyScalar(tempEntr)), vol2.entranceNormal, this);

                        tempEntr = vol2.exit;
                        tempObjEntr = vol2.object;

                        v2Index++;
                    }
                    // old volume ends before vol1 ends
                    else if (tempExit < vol1.exit)
                    {
                        if (tempExit != -FLT_MAX)
                            ray->pushVolumeBoolean(tempEntr, vol2.entrance, tempObjEntr->getNormal(ray->applyScalar(tempEntr)), vol2.entranceNormal, this);

                        tempEntr = vol2.exit;
                        tempExit = vol1.exit;
                        tempObjEntr = vol2.object;
                        tempObjExit = vol1.object;

                        v1Index++; //CHECK
                        v2Index++;
                    }
                    // old volume extends past vol1
                    else
                    {
                        if (tempExit != -FLT_MAX)
                            ray->pushVolumeBoolean(tempEntr, vol2.entrance, tempObjEntr->getNormal(ray->applyScalar(tempEntr)), vol2.entranceNormal, this);

                        tempEntr = vol2.exit;
                        tempObjEntr = vol2.object;

                        v1Index++;  //CHECK
                        v2Index++;
                    }
                }
                // vol1 and vol2 partially overlap
                else if (vol2.exit < vol1.exit)
                {
                    // old volume ends before vol2 begins
                    if (tempExit < vol2.entrance)
                    {
                        if (tempExit != -FLT_MAX)
                            ray->pushVolumeBoolean(tempEntr, tempExit, tempObjEntr->getNormal(ray->applyScalar(tempEntr)), tempObjExit->getNormal(ray->applyScalar(tempExit)), this);

                        tempEntr = vol2.exit;
                        tempExit = vol1.exit;
                        tempObjEntr = vol2.object;
                        tempObjExit = vol1.object;

                        v1Index++;  //CHECK
                        v2Index++;
                    }
                    // old volume ends before vol1 begins
                    else if (tempExit < vol1.entrance)
                    {
                        if (tempExit != -FLT_MAX)
                            ray->pushVolumeBoolean(tempEntr, vol2.entrance, tempObjEntr->getNormal(ray->applyScalar(tempEntr)), vol2.entranceNormal, this);

                        tempEntr = vol2.exit;
                        tempExit = vol1.exit;
                        tempObjEntr = vol2.object;
                        tempObjExit = vol1.object;

                        v1Index++;
                        v2Index++;
                    }
                    // old volume ends before vol2 ends
                    else if (tempExit < vol2.exit)
                    {
                        if (tempExit != -FLT_MAX)
                            ray->pushVolumeBoolean(tempEntr, vol2.entrance, tempObjEntr->getNormal(ray->applyScalar(tempEntr)), vol2.entranceNormal, this);

                        tempEntr = vol2.exit;
                        tempExit = vol1.exit;
                        tempObjEntr = vol2.object;
                        tempObjExit = vol1.object;

                        v1Index++;  //CHECK
                        v2Index++;
                    }
                    // old volume ends before vol1 ends
                    else if (tempExit < vol1.exit)
                    {
                        if (tempExit != -FLT_MAX)
                            ray->pushVolumeBoolean(tempEntr, vol2.entrance, tempObjEntr->getNormal(ray->applyScalar(tempEntr)), vol2.entranceNormal, this);

                        tempEntr = vol2.exit;
                        tempExit = vol1.exit;
                        tempObjEntr = vol2.object;
                        tempObjExit = vol1.object;

                        v1Index++;  //CHECK
                        v2Index++;
                    }
                    // old volume extends past vol1
                    else
                    {
                        if (tempExit != -FLT_MAX)
                            ray->pushVolumeBoolean(tempEntr, vol2.entrance, tempObjEntr->getNormal(ray->applyScalar(tempEntr)), vol2.entranceNormal, this);

                        tempEntr = vol2.exit;
                        tempObjEntr = vol2.object;

                        v1Index++;
                        v2Index++;
                    }
                }
                // vol2 encompasses vol1
                else
                {
                    // old volume ends before vol2 begins
                    if (tempExit < vol2.entrance)
                    {
                        if (tempExit != -FLT_MAX)
                            ray->pushVolumeBoolean(tempEntr, tempExit, tempObjEntr->getNormal(ray->applyScalar(tempEntr)), tempObjExit->getNormal(ray->applyScalar(tempExit)), this);

                        tempEntr = -FLT_MAX;
                        tempExit = -FLT_MAX;
                        tempObjEntr = NULL;
                        tempObjExit = NULL;

                        v1Index++;
                        //v2Index++;    //CHECK
                    }
                    // old volume ends before vol1 begins
                    else if (tempExit < vol1.entrance)
                    {
                        if (tempExit != -FLT_MAX)
                            ray->pushVolumeBoolean(tempEntr, vol2.entrance, tempObjEntr->getNormal(ray->applyScalar(tempEntr)), vol2.entranceNormal, this);

                        tempEntr = -FLT_MAX;
                        tempExit = -FLT_MAX;
                        tempObjEntr = NULL;
                        tempObjExit = NULL;

                        v1Index++;
                        //v2Index++;    //CHECK
                    }
                    // old volume ends before vol1 ends
                    else if (tempExit < vol1.exit)
                    {
                        if (tempExit != -FLT_MAX)
                            ray->pushVolumeBoolean(tempEntr, vol2.entrance, tempObjEntr->getNormal(ray->applyScalar(tempEntr)), vol2.entranceNormal, this);

                        tempEntr = -FLT_MAX;
                        tempExit = -FLT_MAX;
                        tempObjEntr = NULL;
                        tempObjExit = NULL;

                        v1Index++;
                        //v2Index++;    //CHECK
                    }
                    // old volume ends before vol2 ends
                    else if (tempExit < vol2.exit)
                    {
                        if (tempExit != -FLT_MAX)
                            ray->pushVolumeBoolean(tempEntr, vol2.entrance, tempObjEntr->getNormal(ray->applyScalar(tempEntr)), vol2.entranceNormal, this);

                        tempEntr = -FLT_MAX;
                        tempExit = -FLT_MAX;
                        tempObjEntr = NULL;
                        tempObjExit = NULL;

                        v1Index++;
                        //v2Index++;    //CHECK
                    }
                    // old volume extends past vol2
                    else
                    {
                        if (tempExit != -FLT_MAX)
                            ray->pushVolumeBoolean(tempEntr, vol2.entrance, tempObjEntr->getNormal(ray->applyScalar(tempEntr)), vol2.entranceNormal, this);

                        tempEntr = vol2.exit;
                        tempObjEntr = vol2.object;

                        v1Index++;
                        v2Index++;
                    }
                }
            }
        }

        if (tempExit != -FLT_MAX)
            ray->pushVolumeBoolean(tempEntr, tempExit, tempObjEntr->getNormal(ray->applyScalar(tempEntr)), tempObjExit->getNormal(ray->applyScalar(tempExit)), this);
    }
    glm::vec3 getNormal(glm::vec3 intersection)
    {
        glm::vec3 lNormal = leftChild->getNormal(intersection);
        if (lNormal != ZERO_VECTOR)
            return (differenceB) ? -lNormal : lNormal;

        glm::vec3 rNormal = rightChild->getNormal(intersection);
        if (rNormal != ZERO_VECTOR)
            return (differenceB) ? -rNormal : rNormal;

        return ZERO_VECTOR;
    }
    void scale(bool enlarge){}
    void move(glm::vec3 move)
    {
        leftChild->move(move);
        rightChild->move(move);
    }
    void rotate(glm::vec3 rotate){}
    void breakBoolean(std::vector<Object*> *objectVec, int index)
    {
        objectVec->erase(objectVec->begin() + index);
        objectVec->push_back(leftChild);
        objectVec->push_back(rightChild);
        rightChild->differenceB = false;
        delete this;
    }
};
