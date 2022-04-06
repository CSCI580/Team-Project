//
// Created by Peida Han on 04/06/22.
//

#ifndef PROJECT_BOX_H
#define PROJECT_BOX_H

#include "Hittable.h"
#include "Vec3.h"

// 16:51
class Box : public Hittable {
public:
    Box() {}
    Box(const Vec3 &b0, const Vec3 &b1, std::shared_ptr<Material> m)
    { bounds[0] = b0, bounds[1] = b1; };

    virtual bool hit(
            const Ray& r, double t_min, double t_max, hit_record& rec) const override;

public:
    std::shared_ptr<Material> mat_ptr;
    Point3 center;
    double radius;
    Vec3 bounds[2];

};

bool Box::hit(const Ray& r, double t_min, double t_max, hit_record& rec) const {
//    Vec3 oc = r.origin() - center;
//    double a = r.direction().length_squared();
//    double half_b = dot(oc, r.direction());
//    double c = oc.length_squared() - radius*radius;
//
//    double discriminant = half_b*half_b - a*c;
//    if (discriminant < 0) return false;
//    double sqrtd = sqrt(discriminant);
//
//    // Find the nearest root that lies in the acceptable range.
//    double root = (-half_b - sqrtd) / a;
//    if (root < t_min || t_max < root) {
//        root = (-half_b + sqrtd) / a;
//        if (root < t_min || t_max < root)
//            return false;
//    }
//
//    rec.t = root;
//    rec.p = r.at(rec.t);
//    Vec3 outward_normal = (rec.p - center) / radius;
//    rec.set_face_normal(r, outward_normal);
//    rec.mat_ptr = mat_ptr;
    float tmin, tmax, tymin, tymax, tzmin, tzmax;

    tmin = (bounds[r.sign[0]].e[0] - r.orig.e[0]) * r.invdir.e[0];
    tmax = (bounds[1-r.sign[0]].e[0] - r.orig.e[0]) * r.invdir.e[0];
    tymin = (bounds[r.sign[1]].e[1] - r.orig.e[1]) * r.invdir.e[1];
    tymax = (bounds[1-r.sign[1]].e[1] - r.orig.e[1]) * r.invdir.e[1];

    if ((tmin > tymax) || (tymin > tmax))
        return false;

    if (tymin > tmin)
        tmin = tymin;
    if (tymax < tmax)
        tmax = tymax;

    tzmin = (bounds[r.sign[2]].e[2] - r.orig.e[2]) * r.invdir.e[2];
    tzmax = (bounds[1-r.sign[2]].e[2] - r.orig.e[2]) * r.invdir.e[2];

    if ((tmin > tzmax) || (tzmin > tmax))
        return false;

    if (tzmin > tmin)
        tmin = tzmin;
    if (tzmax < tmax)
        tmax = tzmax;

    t = tmin;

    if (t < 0) {
        t = tmax;
        if (t < 0) return false;
    }

    return true;
}

#endif //PROJECT_SPHERE_H
