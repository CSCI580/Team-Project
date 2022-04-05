cmake_minimum_required(VERSION 3.21)
project(Project)

set(CMAKE_CXX_STANDARD 11)

add_executable(Project main.cpp Vec3.h Color.h Ray.h util.h Hittable.h Sphere.h HittableList.h Camera.h)
