
#include "Hittable.h"
#include "Vec3.h"

class translate : public Hittable {
public:
    translate(std::shared_ptr<Hittable> p, const Vec3& displacement)
            : ptr(p), offset(displacement) {}

    virtual bool hit(
            const Ray& r, double t_min, double t_max, hit_record& rec) const override;

    virtual bool bounding_box(float time0, float time1, AABB& output_box) const override;

public:
    std::shared_ptr<Hittable> ptr;
    Vec3 offset;
};

bool translate::hit(const Ray& r, double t_min, double t_max, hit_record& rec) const {
    Ray moved_r(r.origin() - offset, r.direction(), r.time());
    if (!ptr->hit(moved_r, t_min, t_max, rec))
        return false;

    rec.p += offset;
    rec.set_face_normal(moved_r, rec.normal);

    return true;
}

bool translate::bounding_box(float time0, float time1, AABB& output_box) const {
    if (!ptr->bounding_box(time0, time1, output_box))
        return false;

    output_box = AABB(
            output_box.min() + offset,
            output_box.max() + offset);

    return true;
}