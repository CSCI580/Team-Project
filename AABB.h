#include "Hittable.h"
#include "Vec3.h"

inline float ffmin(float a, float b) { return a < b ? a : b; }
inline float ffmax(float a, float b) { return a > b ? a : b; }

// object volumn
class AABB {
public:
    Vec3 _min; //左下角顶点
    Vec3 _max; //右上角顶点
    AABB() {}
    AABB(const Vec3& a, const Vec3& b) { _min = a; _max = b; }

    Vec3 min() const { return _min; }
    Vec3 max() const { return _max; }

    //判断射线是否在t区间
    bool hit(const Ray &r, float tmin, float tmax) const
    {
        for (int a = 0; a < 3; a++)
        {
            float invD = 1.0f / r.direction()[a];
            float t0 = (min()[a] - r.origin()[a]) * invD;
            float t1 = (max()[a] - r.origin()[a]) * invD;
            if (invD < 0.0f)
                std::swap(t0, t1);
            tmax = t1 < tmax ? t1 : tmax;   //F为两者终点的最小值
            tmin = t0 > tmin ? t0 : tmin;   //f为两者起点的最大值
            if (tmax <= tmin)       //F <= f
                return false;
        }
        return true;
    }
};

//计算2个box的包围盒
AABB surrounding_box(AABB box0, AABB box1)
{
    Vec3 small(ffmin(box0.min().x(), box1.min().x()),
               ffmin(box0.min().y(), box1.min().y()),
               ffmin(box0.min().z(), box1.min().z()));
    Vec3 big(ffmax(box0.max().x(), box1.max().x()),
             ffmax(box0.max().y(), box1.max().y()),
             ffmax(box0.max().z(), box1.max().z()));
    return AABB(small, big);
}