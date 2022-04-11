//
// Created by Peida Han on 04/06/22.
//
// reference: https://www.scratchapixel.com/lessons/3d-basic-rendering/minimal-ray-tracer-rendering-simple-shapes/ray-box-intersection
#ifndef PROJECT_BOX_H
#define PROJECT_BOX_H

#include "Hittable.h"
#include "Vec3.h"
#include "AArect.h"

class Box: public Hittable  {
public:
    Box() {}
    Box(const Vec3& p0, const Vec3& p1, std::shared_ptr<Material> m);
    virtual bool hit( const Ray& r, double t0, double t1, hit_record& rec) const;
    virtual bool bounding_box(float t0, float t1, AABB& box) const {
        box =  AABB(pmin, pmax);
        return true; }
    Vec3 pmin, pmax;
    HittableList list;
};

Box::Box(const Vec3& p0, const Vec3& p1, std::shared_ptr<Material> ptr) {
    pmin = p0;
    pmax = p1;
    list.add(std::make_shared<xy_rect>(p0.x(), p1.x(), p0.y(), p1.y(), p1.z(), ptr));
    list.add(std::make_shared<xy_rect>(p0.x(), p1.x(), p0.y(), p1.y(), p0.z(), ptr));

    list.add(std::make_shared<xz_rect>(p0.x(), p1.x(), p0.z(), p1.z(), p1.y(), ptr));
    list.add(std::make_shared<xz_rect>(p0.x(), p1.x(), p0.z(), p1.z(), p0.y(), ptr));

    list.add(std::make_shared<yz_rect>(p0.y(), p1.y(), p0.z(), p1.z(), p1.x(), ptr));
    list.add(std::make_shared<yz_rect>(p0.y(), p1.y(), p0.z(), p1.z(), p0.x(), ptr));

}

bool Box::hit(const Ray& r, double t0, double t1, hit_record& rec) const {
    return list.hit(r, t0, t1, rec);
}

#endif //PROJECT_BOX_H
