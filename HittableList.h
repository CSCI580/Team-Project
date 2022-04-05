//
// Created by Haoyun Zhu on 3/30/22.
//

#ifndef PROJECT_HITTABLELIST_H
#define PROJECT_HITTABLELIST_H

#include "Hittable.h"

#include <memory>
#include <vector>


class HittableList : public Hittable {
public:
    HittableList() {}
    HittableList(std::shared_ptr<Hittable> object) { add(object); }

    void clear() { objects.clear(); }
    void add(std::shared_ptr<Hittable> object) { objects.push_back(object); }

    virtual bool hit(
            const Ray& r, double t_min, double t_max, hit_record& rec) const override;

public:
    std::vector< std::shared_ptr<Hittable> > objects;
};

bool HittableList::hit(const Ray& r, double t_min, double t_max, hit_record& rec) const {
    hit_record temp_rec;
    bool hit_anything = false;
    double closest_so_far = t_max;

    for (const auto& object : objects) {
        if (object->hit(r, t_min, closest_so_far, temp_rec)) {
            hit_anything = true;
            closest_so_far = temp_rec.t;
            rec = temp_rec;
        }
    }

    return hit_anything;
}

#endif //PROJECT_HITTABLELIST_H
