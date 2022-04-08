#include "Hittable.h"
#include "Vec3.h"

class flip_normals : public Hittable {
public:
    flip_normals(Hittable *p) : ptr(p) {}
    virtual bool hit(const Ray& r, double t_min, double t_max, hit_record& rec) const {
        if (ptr->hit(r, t_min, t_max, rec)) {
            rec.normal = -rec.normal;
            return true;
        }
        else
            return false;
    }
    virtual bool bounding_box(float t0, float t1, AABB& box) const {
        return ptr->bounding_box(t0, t1, box);
    }
    Hittable *ptr;
};

class xy_rect : public Hittable {
public:
    xy_rect() {}
    xy_rect(float _x0, float _x1, float _y0, float _y1, float _k, std::shared_ptr<Material> m) : x0(_x0), x1(_x1), y0(_y0), y1(_y1), k(_k), mp(m) {};
    virtual bool hit(const Ray& r, double t_min, double t_max, hit_record& rec) const;

    virtual bool bounding_box(float t0, float t1, AABB& box) const {
        box = AABB(Vec3(x0, y0, k - 0.0001), Vec3(x1, y1, k + 0.0001));
        return true;
    }
    std::shared_ptr<Material> mp;
    float x0, x1, y0, y1, k;
};

bool xy_rect::hit(const Ray& r, double t0, double t1, hit_record& rec) const {
    float t = (k - r.origin().z()) / r.direction().z();
    if (t < t0 || t > t1)        return false;
    float x = r.origin().x() + t * r.direction().x();
    float y = r.origin().y() + t * r.direction().y();
    if (x < x0 || x > x1 || y < y0 || y > y1)         return false;
    rec.u = (x - x0) / (x1 - x0);
    rec.v = (y - y0) / (y1 - y0);
    rec.t = t;
    rec.mat_ptr = mp;
    rec.p = r.at(t);
    rec.normal = Vec3(0, 0, 1);	//默认法线为(0,0,1)，但会有潜在的朝向问题
    return true;
}



class xz_rect : public Hittable {
public:
    xz_rect() {}
    xz_rect(float _x0, float _x1, float _z0, float _z1, float _k, std::shared_ptr<Material> m) :x0(_x0), x1(_x1), z0(_z0), z1(_z1), k(_k), mp(m) {};
    virtual bool hit(const Ray& r, double t_min, double t_max, hit_record& rec) const{
        auto t = (k - r.origin().y()) / r.direction().y();
        if (t < t_min || t > t_max)
            return false;
        auto x = r.origin().x() + t * r.direction().x();
        auto z = r.origin().z() + t * r.direction().z();
        if (x < x0 || x > x1 || z < z0 || z > z1)
            return false;
        rec.u = (x - x0) / (x1 - x0);
        rec.v = (z - z0) / (z1 - z0);
        rec.t = t;
        rec.normal = Vec3(0, 1, 0);
        rec.mat_ptr = mp;
        rec.p = r.at(t);
        return true;
    };

    virtual bool bounding_box(float time0, float time1, AABB& output_box) const override {
        // The bounding box must have non-zero width in each dimension, so pad the Y
        // dimension a small amount.
        output_box = AABB(Vec3(x0, k - 0.0001, z0), Vec3(x1, k + 0.0001, z1));
        return true;
    }
    std::shared_ptr<Material> mp;
    double x0, x1, z0, z1, k;
};


class yz_rect : public Hittable {
public:

    yz_rect(float _y0, float _y1, float _z0, float _z1, float _k,std::shared_ptr<Material> m): y0(_y0), y1(_y1), z0(_z0), z1(_z1), k(_k), mp(m) {};

    virtual bool hit(const Ray& r, double t_min, double t_max, hit_record& rec) const {
        auto t = (k - r.origin().x()) / r.direction().x();
        if (t < t_min || t > t_max)
            return false;
        auto y = r.origin().y() + t * r.direction().y();
        auto z = r.origin().z() + t * r.direction().z();
        if (y < y0 || y > y1 || z < z0 || z > z1)
            return false;
        rec.u = (y - y0) / (y1 - y0);
        rec.v = (z - z0) / (z1 - z0);
        rec.t = t;
        rec.normal = Vec3(1, 0, 0);

        rec.mat_ptr = mp;
        rec.p = r.at(t);
        return true;
    };

    virtual bool bounding_box(float time0, float time1, AABB& output_box) const override {
        // The bounding box must have non-zero width in each dimension, so pad the X
        // dimension a small amount.
        output_box = AABB(Vec3(k - 0.0001, y0, z0), Vec3(k + 0.0001, y1, z1));
        return true;
    }

    std::shared_ptr<Material> mp;
    double y0, y1, z0, z1, k;
};
