//
// Created by Haoyun Zhu on 3/30/22.
//

#ifndef PROJECT_VEC3_H
#define PROJECT_VEC3_H

#include <cmath>
#include <iostream>
#include "util.h"

#define Point3 Vec3
#define Color Vec3

class Vec3 {
public:
    Vec3() {
        e[0] = 0; e[1] = 0; e[2] = 0;
    }
    Vec3(double e0, double e1, double e2){
        e[0] = e0; e[1] = e1; e[2] = e2;
    }

    double x() const { return e[0]; }
    double y() const { return e[1]; }
    double z() const { return e[2]; }

    Vec3 operator-() const { return Vec3(-e[0], -e[1], -e[2]); }
    double operator[](int i) const { return e[i]; }
    double& operator[](int i) { return e[i]; }

    Vec3& operator+=(const Vec3 &v) {
        e[0] += v.e[0];
        e[1] += v.e[1];
        e[2] += v.e[2];
        return *this;
    }

    Vec3& operator*=(const double t) {
        e[0] *= t;
        e[1] *= t;
        e[2] *= t;
        return *this;
    }

    Vec3& operator/=(const double t) {
        return *this *= 1/t;
    }

    double length() const {
        return sqrt(length_squared());
    }

    double length_squared() const {
        return e[0]*e[0] + e[1]*e[1] + e[2]*e[2];
    }

    inline static Vec3 random() {
        return Vec3(random_double(), random_double(), random_double());
    }

    inline static Vec3 random(double min, double max) {
        return Vec3(random_double(min,max), random_double(min,max), random_double(min,max));
    }

    inline bool near_zero() const {
        // Return true if the vector is close to zero in all dimensions.
        const double s = EPSILON;
        return (fabs(e[0]) < s) && (fabs(e[1]) < s) && (fabs(e[2]) < s);
    }

public:
    double e[3];
};
// vec3 Utility Functions

inline std::ostream& operator<<(std::ostream &out, const Vec3 &v) {
    return out << v.e[0] << ' ' << v.e[1] << ' ' << v.e[2];
}

inline Vec3 operator+(const Vec3 &u, const Vec3 &v) {
    return Vec3(u.e[0] + v.e[0], u.e[1] + v.e[1], u.e[2] + v.e[2]);
}

inline Vec3 operator-(const Vec3 &u, const Vec3 &v) {
    return Vec3(u.e[0] - v.e[0], u.e[1] - v.e[1], u.e[2] - v.e[2]);
}

inline Vec3 operator*(const Vec3 &u, const Vec3 &v) {
    return Vec3(u.e[0] * v.e[0], u.e[1] * v.e[1], u.e[2] * v.e[2]);
}

inline Vec3 operator*(double t, const Vec3 &v) {
    return Vec3(t*v.e[0], t*v.e[1], t*v.e[2]);
}

inline Vec3 operator*(const Vec3 &v, double t) {
    return t * v;
}

inline Vec3 operator/(Vec3 v, double t) {
    return (1/t) * v;
}

inline double dot(const Vec3 &u, const Vec3 &v) {
    return u.e[0] * v.e[0]
           + u.e[1] * v.e[1]
           + u.e[2] * v.e[2];
}

inline Vec3 cross(const Vec3 &u, const Vec3 &v) {
    return Vec3(u.e[1] * v.e[2] - u.e[2] * v.e[1],
                u.e[2] * v.e[0] - u.e[0] * v.e[2],
                u.e[0] * v.e[1] - u.e[1] * v.e[0]);
}

inline Vec3 normalized(Vec3 v) {
    double l = v.length();
    if (equals(l,0)) return Vec3(0,0,0);
    else return v / v.length();
}

inline Vec3 unit_vector(Vec3 v) {
    double l = v.length();
    if (equals(l,0)) return Vec3(0,0,0);
    else return v / v.length();
}

inline Vec3 random_unit_vector() {
    double theta = 2 * M_PI * random_double(0,1);
    double phi = acos(1 - 2 * random_double(0,1));
    double x = sin(phi) * cos(theta);
    double y = sin(phi) * sin(theta);
    double z = cos(phi);
    return Vec3(x,y,z);
}

inline Vec3 random_in_unit_sphere() {
    while (true) {
        Vec3 p = Vec3::random(-1,1);
        if (p.length_squared() >= 1) continue;
        return p;
    }
}

Vec3 random_in_unit_disk() {
    while (true) {
        auto p = Vec3(random_double(-1,1), random_double(-1,1), 0);
        if (p.length_squared() >= 1) continue;
        return p;
    }
}

Vec3 reflect(const Vec3& v, const Vec3& n) {
    return v - 2*dot(v,n)*n;
}

Vec3 refract(const Vec3& uv, const Vec3& n, double etai_over_etat) {
    double cos_theta = fmin(dot(-uv, n), 1.0);
    Vec3 r_out_perp =  etai_over_etat * (uv + cos_theta*n);
    Vec3 r_out_parallel = -sqrt(fabs(1.0 - r_out_perp.length_squared())) * n;
    return r_out_perp + r_out_parallel;
}

#endif //PROJECT_VEC3_H
