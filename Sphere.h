//
// Created by Haoyun Zhu on 3/30/22.
//

#ifndef PROJECT_SPHERE_H
#define PROJECT_SPHERE_H

#include "Hittable.h"
#include "Vec3.h"

// 16:51
class Sphere : public Hittable {
public:
    Sphere() {}
    Sphere(Point3 cen, double r, std::shared_ptr<Material> m)
    : center(cen), radius(r), mat_ptr(m) {};

    virtual bool hit(
            const Ray& r, double t_min, double t_max, hit_record& rec) const override;
    virtual bool bounding_box(float t0, float t1, AABB& box) const override;
public:
    std::shared_ptr<Material> mat_ptr;
    Point3 center;
    double radius;
};

bool Sphere::hit(const Ray& r, double t_min, double t_max, hit_record& rec) const {
    Vec3 oc = r.origin() - center;
    double a = r.direction().length_squared();
    double half_b = dot(oc, r.direction());
    double c = oc.length_squared() - radius*radius;

    double discriminant = half_b*half_b - a*c;
    if (discriminant < 0) return false;
    double sqrtd = sqrt(discriminant);

    // Find the nearest root that lies in the acceptable range.
    double root = (-half_b - sqrtd) / a;
    if (root < t_min || t_max < root) {
        root = (-half_b + sqrtd) / a;
        if (root < t_min || t_max < root)
            return false;
    }

    rec.t = root;
    rec.p = r.at(rec.t);
    Vec3 outward_normal = (rec.p - center) / radius;
    rec.set_face_normal(r, outward_normal);
    rec.mat_ptr = mat_ptr;

    return true;
}
bool Sphere::bounding_box(float t0, float t1, AABB& box) const {
    box = AABB(
            center - Vec3(radius, radius, radius),
            center + Vec3(radius, radius, radius));
    return true;
}
#endif //PROJECT_SPHERE_H
