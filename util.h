//
// Created by Haoyun Zhu on 3/30/22.
//

#ifndef PROJECT_UTIL_H
#define PROJECT_UTIL_H

#include <cmath>
#include <memory>
#include <random>

// Constants
const double EPSILON = 1e-10;
const double PI = 3.1415926535897932385;
const double INF = std::numeric_limits<double>::infinity();


// Utility Functions
inline bool equals (double x, double y) {
    return fabs(x-y) < EPSILON;
}

inline double degrees_to_radians(double degrees) {
    return degrees * PI / 180.0;
}

inline double random_double() {
    static std::uniform_real_distribution<double> distribution(0.0, 1.0);
    static std::mt19937 generator;
    return distribution(generator);
}

inline double random_double(double min, double max) {
    static std::uniform_real_distribution<double> distribution(min, max);
    static std::mt19937 generator;
    return distribution(generator);
}

inline double clamp(double x, double min, double max) {
    if (x < min) return min;
    if (x > max) return max;
    return x;
}

// Common Headers
#include "Ray.h"
#include "Vec3.h"

#endif //PROJECT_UTIL_H
