cmake_minimum_required(VERSION 3.21)
project(Project)

set(CMAKE_CXX_STANDARD 11)

add_executable(openMP main.cpp)
FIND_PACKAGE(OpenMP REQUIRED)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()


## added -fopenmp
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fopenmp")
#
#set(SOURCE_FILES main.cpp)
#add_executable(openmp_test_clion ${SOURCE_FILES})


add_executable(Project main.cpp Vec3.h Color.h Ray.h util.h Hittable.h Sphere.h HittableList.h Camera.h)
