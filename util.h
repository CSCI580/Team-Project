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
    // Returns a random real in [0,1).
    static std::random_device rd;
    static std::mt19937 generator(20000905 << 2); //Mersenne Twister 19937 generator
    static std::uniform_real_distribution<double> distribution(0.0, 1.0);
    return distribution(generator);
}

inline double random_double(double min, double max) {
    return min + (max - min) * random_double();
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
