#pragma once
#include <algorithm>
#include <glm/glm.hpp>
#include <vector>
#include <utility>
#include <glm/gtx/rotate_vector.hpp>
#include <complex>

// Mathematical values
#define PI			3.14159265359f
#define identity	glm::mat4(1.f)
#define FLOAT_ERROR 1e-4
#define radToDeg    (PI / 180)

#define BASE_CENTER glm::vec3(0.f);
#define BASE_RADIUS .5f
#define BASE_PHONG  50
#define BASE_COLOUR glm::vec3(0.f, 1.f, 0.f)

#define SCALE_CHANGE .01f
#define MIN_SCALE .1f

#define XAXIS glm::vec3(1.f, 0.f, 0.f)
#define YAXIS glm::vec3(0.f, 1.f, 0.f)
#define ZAXIS glm::vec3(0.f, 0.f, 1.f)

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
    std::vector<Volume> getVolume() { return volumes; }
    void pushVolume(float s1, float s2, Object *obj)
    {
        volumes.push_back(Volume(s1, s2, obj));
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
    Ray generateTempRay(Ray *ray, glm::vec3 rotation)
    {
		glm::vec3 originRotated = glm::rotateX(glm::rotateY(glm::rotateZ(ray->origin, -rotation.z), -rotation.y), -rotation.x);
        glm::vec3 dirRotated = glm::rotateX(glm::rotateY(glm::rotateZ(ray->direction, -rotation.z), -rotation.y), -rotation.x);
        return Ray(originRotated - center, dirRotated);
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

// objects only used by other objects
struct Plane : Object
{
    glm::vec3 normal;
    Plane(glm::vec3 point, glm::vec3 normal) :
        normal(normalize(normal))
    {
        center = point;
        colour = BASE_COLOUR;
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

// the actual objects
struct Sphere : Object
{
    Sphere()
    {
        radius = BASE_RADIUS;
        center = BASE_CENTER;
        colour = BASE_COLOUR;
        phong = BASE_PHONG;
    }

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
    // variables
    Plane* planes[6];
    glm::vec3 top, side, front;
    float radi;

    // constructors
    Cube()
    {
        radius = BASE_RADIUS;
        center = BASE_CENTER;
        top = YAXIS;
        side = XAXIS;
        front = ZAXIS;
        colour = BASE_COLOUR;
        phong = BASE_PHONG;

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
    glm::vec3 rotation;
    float R, r;

    Torus()
    {
        center = BASE_CENTER;
        radius = BASE_RADIUS;
        colour = BASE_COLOUR;
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
        glm::vec3 originRotated = glm::rotateX(glm::rotateY(glm::rotateZ(ray->origin, -rotation.z), -rotation.y), -rotation.x);
        glm::vec3 dirRotated = glm::rotateX(glm::rotateY(glm::rotateZ(ray->direction, -rotation.z), -rotation.y), -rotation.x);
        Ray tempRay(originRotated - center, dirRotated);

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
            {
               scalars.push_back(roots[i]);
            }
        }




        std::sort(scalars.begin(), scalars.end());
        switch (scalars.size())
        {
        case(2):
            ray->volumes.push_back(Volume(scalars[0], scalars[1], this));
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
    void scale(bool enlarge) 
    {
        R += (enlarge) ? SCALE_CHANGE : -SCALE_CHANGE;
        r = R * .25f;
    }
    void move(glm::vec3 move) { center += move; }
    void rotate(glm::vec3 rotate) { rotation += rotate; }
    void select() { selected = true; }
    void deselect() { selected = false; }
    void breakBoolean(std::vector<Object*> *objectVec, int index) {}
};

struct Cylinder : Object
{
    #define BASE_LENGTH .5f
    float length;
    glm::vec3 rotation = glm::vec3(0.f);
    Plane *topPlane, *bottomPlane;

    Cylinder()
    {
        radius = BASE_RADIUS;
        center = BASE_CENTER;
        colour = BASE_COLOUR;
        phong = BASE_PHONG;
        length = BASE_LENGTH;

        topPlane = new Plane(center + (length * YAXIS), YAXIS);
        bottomPlane = new Plane(center - (length * YAXIS), -YAXIS);
    }
    ~Cylinder() { delete topPlane; delete bottomPlane; }

    void getVolume(Ray *ray)
    {
        Ray tempRay = generateTempRay(ray, rotation);

        std::vector<float> scalars;

        // get and check intersections with planes
        float tempScalar = topPlane->getIntersection(&tempRay);
        if (tempScalar != NULL)
        {
            glm::vec3 intersection = tempRay.applyScalar(tempScalar);
            if (abs(distance(intersection, glm::vec3(center.x, intersection.y, center.z)) <= radius))
                scalars.push_back(tempScalar);
        }
        tempScalar = bottomPlane->getIntersection(&tempRay);
        if (tempScalar != NULL)
        {
            glm::vec3 intersection = tempRay.applyScalar(tempScalar);
            if (abs(distance(intersection, glm::vec3(center.x, intersection.y, center.z)) <= radius))
                scalars.push_back(tempScalar);
        }

        // if the ray collides perfectly with both planes then it is inside the cylinder, and we dont need to check the cylinder
        if (scalars.size() == 2)
        {
            ray->volumes.push_back(Volume(glm::min(scalars[0], scalars[1]), glm::max(scalars[0], scalars[1]), this));
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
        ray->volumes.push_back(Volume(scalars.front(), scalars.back(), this));
    }
    glm::vec3 getNormal(glm::vec3 intersection)
    {
		glm::vec2 inter = glm::vec2(intersection.x, intersection.z);
		glm::vec3 norm(0.f);
			
		if (abs(dot(normalize(intersection - topPlane->center), topPlane->normal)) <= FLOAT_ERROR)
            norm = topPlane->normal;
            
		else if (abs(dot(normalize(intersection - bottomPlane->center), bottomPlane->normal)) <= FLOAT_ERROR)
            norm = bottomPlane->normal;
           
		else if (distance(inter, glm::vec2(center.x, center.z)) - radius <= FLOAT_ERROR)
			norm = normalize(intersection - glm::vec3(center.x, intersection.y, center.z));
			
		else 
			return glm::vec3(0.f);
			
        return norm;

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
    void select() 
    { 
        selected = true; 
        topPlane->select();
        bottomPlane->select();
    }
    void deselect() 
    { 
        selected = false; 
        topPlane->deselect();
        bottomPlane->deselect();
    }
    void breakBoolean(std::vector<Object*> *objectVec, int index) {}
};

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// the combinations
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
        float v1Index = 0;
        float v2Index = 0;

        // temporary storge for enterences and exits
        // initialize the floats to smallest possible
        float tempEntr = -FLT_MAX;
        float tempExit = -FLT_MAX;

        while (v1Index != v1.size() && v2Index != v2.size())
        {
            // define as new volume for easier access
            Volume vol1 = v1[v1Index], vol2 = v2[v2Index];

            // volume 1 is closer to the volume
            if (vol1.entrance < vol2.entrance)
            {
                if (tempExit < vol1.entrance)
                {
                    // old volume is not involved. push it to the ray
                    if (tempExit != -FLT_MAX)
                        ray->pushVolume(tempEntr, tempExit, this);

                    // volume 1 is not involved in the intersection
                    if (vol1.exit < vol2.entrance)
                    {
                        // reset tempEntr and tempExit to defaults
                        tempEntr = -FLT_MAX;
                        tempExit = -FLT_MAX;
                        // toss v1, and restart volume search
                        v1Index++;
                        continue;
                    }
                    // volume 2 entrance is contained within volume 1. reassign th volume
                    else
                    {
                        tempEntr = vol2.entrance;
                        tempExit = glm::min(vol1.exit, vol2.exit);
                        v1Index++;
                        v2Index++;
                    }
                }
                // volume1s entrance is contained in the volume
                else
                {
                    v1Index++;// v1Index is done

                    if (tempExit < vol2.entrance)
                    {
                        if (tempExit != -FLT_MAX)
                            ray->pushVolume(tempEntr, tempExit, this);

                        // reset volume to defaults
                        tempEntr = -FLT_MAX;
                        tempExit = -FLT_MAX;
                    }
                }
            }
            // volume 2 is closer to the volume
            else
            {
                if (tempExit < vol2.entrance)
                {
                    // old volume is not involved. push it to the ray
                    if (tempExit != -FLT_MAX)
                        ray->pushVolume(tempEntr, tempExit, this);

                    // volume 1 is not involved in the intersection
                    if (vol2.exit < vol1.entrance)
                    {
                        // reset tempEntr and tempExit to defaults
                        tempEntr = -FLT_MAX;
                        tempExit = -FLT_MAX;
                        // toss v1, and restart volume search
                        v2Index++;
                        continue;
                    }
                    // volume 2 entrance is contained within volume 1. reassign th volume
                    else
                    {
                        tempEntr = vol1.entrance;
                        tempExit = glm::min(vol2.exit, vol1.exit);
                        v1Index++;
                        v2Index++;
                    }
                }
                // volume2s entrance is contined in the volume
                else
                {
                    v2Index++;// v1Index is done

                    if (tempExit < vol1.entrance)
                    {
                        if (tempExit != -FLT_MAX)
                            ray->pushVolume(tempEntr, tempExit, this);

                        // reset volume to defaults
                        tempEntr = -FLT_MAX;
                        tempExit = -FLT_MAX;
                    }
                }
            }
        }

        // push the remaining volume
        if (tempExit != -FLT_MAX)
            ray->pushVolume(tempEntr, tempExit, this);




        /*
        Ray tempRay(ray->origin, ray->direction);
        leftObject->getVolume(&tempRay);

        bool obj1Torus = false; // need to know order of objects if one i a torus
        if (tempRay.volumes.size() == 2) obj1Torus = true;

        rightObject->getVolume(&tempRay);

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
        */
    }
    glm::vec3 getNormal(glm::vec3 intersection)
    {
        if (distance(intersection, leftChild->center) - leftChild->radius <= FLOAT_ERROR)
            return leftChild->getNormal(intersection);
        else
            return rightChild->getNormal(intersection);
    }
    void scale(bool enlarge)
    {
        leftChild->scale(enlarge);
        rightChild->scale(enlarge);
    }
    void move(glm::vec3 move)
    {
        leftChild->move(move);
        rightChild->move(move);
    }
    void rotate(glm::vec3 rotate)
    {
        leftChild->rotate(rotate);
        rightChild->rotate(rotate);
    }
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

        // handle the cases where one of the lists are empty
        if (v1.size() == 0)
        {
            for (Volume vol : v2)
                ray->pushVolume(vol.entrance, vol.exit, this);
            return;
        }
        if (v2.size() == 0)
        {
            for (Volume vol : v1)
                ray->pushVolume(vol.entrance, vol.exit, this);
            return;
        }


        // intex trackers for each list
        float v1Index = 0;
        float v2Index = 0;

        // temporary storge for enterences and exits
        // initialize the floats to smallest possible
        float tempEntr = -FLT_MAX;
        float tempExit = -FLT_MAX;

        while (v1Index != v1.size() && v2Index != v2.size())
        {
            // define as new volume for easier access
            Volume vol1 = v1[v1Index], vol2 = v2[v2Index];

            // volume 1 is closer than volume 2
            if (vol1.entrance < vol2.entrance)
            { 
                // current volume ends before new volume starts
                if (tempExit < vol1.entrance)
                {
                    if (tempExit != -FLT_MAX)
                        ray->pushVolume(tempEntr, tempExit, this);
                    tempEntr = vol1.entrance;
                    tempExit = vol1.exit;
                }
                // old volume starts before vol1 start and ends before vol1 ends. Union so add on vol 2
                else if (tempExit < vol1.exit)
                    tempExit = vol1.exit;

                // old volume ends before volume 2 starts
                if (tempExit < vol2.entrance)
                {
                    ray->pushVolume(tempEntr, tempExit, this);
                    tempEntr = vol2.entrance;
                    tempExit = vol2.exit;
                }
                // old volume starts before vol2 start and ends before vol2 ends. Union so add on vol 2
                else if (tempExit < vol2.exit)
                    tempExit = vol2.exit;
            }
            // volume 2 is closer than volume 1
            else
            {
                // current volume ends before new volume starts
                if (tempExit < vol2.entrance)
                {
                    if (tempExit != -FLT_MAX)
                        ray->pushVolume(tempEntr, tempExit, this);
                    tempEntr = vol2.entrance;
                    tempExit = vol2.exit;
                }
                // old volume starts before vol1 start and ends before vol1 ends. Union so add on vol 2
                else if (tempExit < vol2.exit)
                    tempExit = vol2.exit;

                // old volume ends before volume 2 starts
                if (tempExit < vol1.entrance)
                {
                    ray->pushVolume(tempEntr, tempExit, this);
                    tempEntr = vol1.entrance;
                    tempExit = vol1.exit;
                }
                // old volume starts before vol2 start and ends before vol2 ends. Union so add on vol 2
                else if (tempExit < vol1.exit)
                    tempExit = vol1.exit;
            }

            // we havce checked both v1Index and v2Index. increment both
            v1Index++;
            v2Index++;
        }

        // push the remaining volume
        if (tempExit != -FLT_MAX)
            ray->pushVolume(tempEntr, tempExit, this);

        // if at the end of the while loop, either list has volumes remaining, since this is union, add these volumes to the ray
        for (v1Index; v1Index < v1.size(); v1Index++)
            ray->pushVolume(v1[v1Index].entrance, v1[v1Index].exit, this);
        for (v2Index; v2Index < v2.size(); v2Index++)
            ray->pushVolume(v2[v2Index].entrance, v2[v2Index].exit, this);
    }
    glm::vec3 getNormal(glm::vec3 intersection)
    {
        if (leftChild->getNormal(intersection) != glm::vec3(0.f))
            return leftChild->getNormal(intersection);
        else
            return rightChild->getNormal(intersection);
    }
    void scale(bool enlarge)
    {
        leftChild->scale(enlarge);
        rightChild->scale(enlarge);
    }
    void move(glm::vec3 move)
    {
        leftChild->move(move);
        rightChild->move(move);
    }
    void rotate(glm::vec3 rotate)
    {
        leftChild->rotate(rotate);
        rightChild->rotate(rotate);
    }
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

struct Difference : Object
{
    Object *leftChild, *rightChild;
    Difference(Object *obj1, Object *obj2) : leftChild(obj1), rightChild(obj2)
    {
        phong = (obj1->phong + obj2->phong) / 2.f;
        colour = (obj1->colour + obj2->colour) / 2.f;
        center = (obj1->center + obj2->center) / 2.f;
        selected = false;
    }

    void getVolume(Ray *ray)
    {
        Ray tempRay(ray->origin, ray->direction);



        leftChild->getVolume(&tempRay);
        // if there is no collision with A in A - B then there is no collision with the difference
        if (tempRay.volumes.size() == 0) return;

        // to toggle special conditions for torus
        bool obj1Torus = false;
        if (tempRay.volumes.size() == 2) obj1Torus = true;


        rightChild->getVolume(&tempRay);


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
        if (abs(distance(intersection, center) - leftChild->radius) <= FLOAT_ERROR)
            return leftChild->getNormal(intersection);
        else
            return -rightChild->getNormal(intersection);
    }
    void scale(bool enlarge)
    {
        leftChild->scale(enlarge);
        rightChild->scale(enlarge);
    }
    void move(glm::vec3 move)
    {
        leftChild->move(move);
        rightChild->move(move);
    }
    void rotate(glm::vec3 rotate)
    {
        leftChild->rotate(rotate);
        rightChild->rotate(rotate);
    }
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
