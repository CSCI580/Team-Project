/*//
// Created by Haoyun Zhu on 3/30/22.
//

#ifndef PROJECT_HITTABLELIST_H
#define PROJECT_HITTABLELIST_H

#include "Hittable.h"

#include <memory>
#include <vector>

class HittableList : public Hittable
{
public:
    Hittable **list;
    int list_size;

    HittableList() {}
   *//* HittableList(Hittable **l, int n)
    {
        list = l;
        list_size = n;
    }*//*
    HittableList(std::shared_ptr<Hittable> object) {
        add(object);
    }
    void add(std::shared_ptr<Hittable> object) {
        objects.push_back(object);
    }

    //如果命中了，命中记录保存到rec
    virtual bool hit(const Ray &r, double t_min, double t_max, hit_record &rec) const
    {
        hit_record temp_rec;
        bool hit_anything = false;
        double closest_so_far = t_max; //记录目前最近的t值
        for (int i = 0; i < list_size; i++)
        {
            if (list[i]->hit(r, t_min, closest_so_far, temp_rec))
            {
                hit_anything = true;
                closest_so_far = temp_rec.t;
                rec = temp_rec; //只记录打到的最近的球
            }
        }
        return hit_anything;
    }

    bool bounding_box(float t0, float t1, AABB &box) const
    {
        if (list_size < 1)
            return false;
        AABB temp_box;
        bool first_true = list[0]->bounding_box(t0, t1, temp_box);
        if (!first_true)
            return false;
        else
            box = temp_box;
        for (int i = 1; i < list_size; i++)
        {
            if (list[i]->bounding_box(t0, t1, temp_box))
            {
                box = surrounding_box(box, temp_box);
            }
            else
                return false;
        }
        return true;
    }
};

#endif //PROJECT_HITTABLELIST_H*/
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
    virtual bool bounding_box(float t0, float t1, AABB& box) const;
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
bool HittableList::bounding_box(float t0, float t1, AABB &box) const
{
    if (objects.size() < 1)
        return false;
    AABB temp_box;
    bool first_true = objects[0]->bounding_box(t0, t1, temp_box);
    if (!first_true)
        return false;
    else
        box = temp_box;
    for (int i = 1; i < objects.size(); i++)
    {
        if (objects[i]->bounding_box(t0, t1, temp_box))
        {
            box = surrounding_box(box, temp_box);
        }
        else
            return false;
    }
    return true;
}
#endif //PROJECT_HITTABLELIST_H
