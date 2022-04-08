//
// Created by Haoyun Zhu on 3/30/22.
//

#ifndef PROJECT_RAY_H
#define PROJECT_RAY_H
#include "Vec3.h"


class Ray {
public:
    Ray() {}
    Ray(const Point3& origin, const Vec3& direction)
            : orig(origin), dir(direction)
    {
        invdir = direction;
        invdir.e[0] = 1 /invdir.e[0];
        invdir.e[1] = 1 /invdir.e[1];
        invdir.e[2] = 1 /invdir.e[2];

        sign[0] = (invdir.e[0] < 0);
        sign[1] = (invdir.e[1] < 0);
        sign[2] = (invdir.e[2] < 0);
    }

    Point3 origin() const  { return orig; }
    Vec3 direction() const { return dir; }

    Point3 at(double t) const {
        return orig + t*dir;
    }

public:
    Point3 orig;
    Vec3 dir;
    Vec3 invdir;
    int sign[3];

};

#endif //PROJECT_RAY_H
