
#include "Hittable.h"
#include "Vec3.h"

class translate : public Hittable {
public:
    translate(Hittable *p, const Vec3& displacement)
            : ptr(p), offset(displacement) {}
    virtual bool hit(
            const Ray& r, double t_min, double t_max, hit_record& rec) const;
    virtual bool bounding_box(float t0, float t1, AABB& box) const;
    Hittable *ptr;
    Vec3 offset;
};

bool translate::hit(const Ray& r, double t_min, double t_max, hit_record& rec) const {
    Ray moved_r(r.origin() - offset, r.direction());
    if (ptr->hit(moved_r, t_min, t_max, rec)) {
        rec.p += offset;
        return true;
    }
    else
        return false;
}

bool translate::bounding_box(float t0, float t1, AABB& box) const {
    if (ptr->bounding_box(t0, t1, box)) {
        box = AABB(box.min() + offset, box.max() + offset);
        return true;
    }
    else
        return false;
}
