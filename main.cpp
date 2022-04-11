#include "util.h"
#include "Color.h"
#include "Hittable.h"
#include "HittableList.h"
#include "Sphere.h"
#include "Camera.h"
#include "Box.h"
#include "Translate.h"
#include "Rotate_Y.h"
#include "BVH_Node.h"

#include <chrono>
#include <memory>
#include <iostream>
#include <fstream>

#define MULTITHREAD 1
#if MULTITHREAD
#include <omp.h>
#endif

HittableList cornell_box() {
    HittableList world;
    auto red = std::make_shared<lambertian>(Vec3(0.65, 0.05, 0.05));
    auto white = std::make_shared<lambertian>(Vec3(0.73, 0.73, 0.73));
    auto green = std::make_shared<lambertian>(Vec3(0.12, 0.45, 0.15));
    auto clear_mat = std::make_shared<dielectric>(1.4);

    world.add(std::make_shared<flip_normals>(new yz_rect(0, 555, 0, 555, 555, green)));
    world.add(std::make_shared<yz_rect>(0, 555, 0, 555, 0, red));
    world.add(std::make_shared<xz_rect>(213, 343, 227, 332, 554, white));
    world.add(std::make_shared<flip_normals>(new xz_rect(0, 555, 0, 555, 555, white)));
    world.add(std::make_shared<xz_rect>(0, 555, 0, 555, 0, white));
    world.add(std::make_shared<flip_normals>(new xy_rect(0, 555, 0, 555, 555, white)));
    // long box
    std::shared_ptr<Hittable> box1 = std::make_shared<Box>(Vec3(0, 0, 0), Vec3(165, 330, 165), white);
    box1 = std::make_shared<rotate_y>(box1, 15);
    box1 = std::make_shared<translate>(box1, Vec3(265,0,295));
    // small box
    std::shared_ptr<Hittable> box2 = std::make_shared<Box>(Vec3(0,0,0), Vec3(165,165,165), clear_mat);
    box2 = std::make_shared<rotate_y>(box2, -18);
    box2 = std::make_shared<translate>(box2, Vec3(130,0,65));

    world.add(box1);
    world.add(box2);

    return world;
}
HittableList random_scene() {
    HittableList world;
    auto white_mat = std::make_shared<lambertian>(Vec3(.73, .73, .73));
    auto clear_mat = std::make_shared<dielectric>(1.4);
    auto metal_mat = std::make_shared<metal>(Color(0.7, 0.6, 0.5), 0.0);
    auto ground_material = std::make_shared<lambertian>(Color(0.5, 0.5, 0.5));
    HittableList boxes1;
//    // long box
//    std::shared_ptr<Hittable> box1 = std::make_shared<Box>(Vec3(0, 0, 0), Vec3(1, 2, 1), white_mat);
//    box1 = std::make_shared<rotate_y>(box1, 15);
//    box1 = std::make_shared<translate>(box1, Vec3(-4, 0, 0));
    // small box
    std::shared_ptr<Hittable> box2 = std::make_shared<Box>(Vec3(0,0,0), Vec3(1,1,1), clear_mat);
//    box2 = std::make_shared<rotate_y>(box2, 0);
    box2 = std::make_shared<translate>(box2, Vec3(0, 0, 0));

//    world.add(box1);
    world.add(box2);

    world.add(std::make_shared<Sphere>(Point3(0,-1000,0), 1000, ground_material));

    for (int a = -11; a < 11; a++) {
        for (int b = -11; b < 11; b++) {
            double choose_mat = random_double();
            Point3 center(a + 0.9*random_double(), 0.2, b + 0.9*random_double());

            if ((center - Point3(4, 0.2, 0)).length() > 0.9) {
                std::shared_ptr<Material> sphere_material;

                if (choose_mat < 0.8) {
                    // diffuse
                    auto albedo = Color::random() * Color::random();
                    sphere_material = std::make_shared<lambertian>(albedo);
                    world.add(std::make_shared<Sphere>(center, 0.2, sphere_material));
                } else if (choose_mat < 0.95) {
                    // metal
                    auto albedo = Color::random(0.5, 1);
                    auto fuzz = random_double(0, 0.5);
                    sphere_material = std::make_shared<metal>(albedo, fuzz);
                    world.add(std::make_shared<Sphere>(center, 0.2, sphere_material));
                } else {
                    // glass
                    sphere_material = std::make_shared<dielectric>(1.5);
                    world.add(std::make_shared<Sphere>(center, 0.2, sphere_material));
                }
            }
        }
    }

 /*   auto material1 = std::make_shared<dielectric>(1.5);
    world.add(std::make_shared<Sphere>(Point3(0, 1, 0), 1.0, material1));

    auto material2 = std::make_shared<lambertian>(Color(0.4, 0.2, 0.1));
    world.add(std::make_shared<Sphere>(Point3(-4, 1, 0), 1.0, material2));

    auto material3 = std::make_shared<metal>(Color(0.7, 0.6, 0.5), 0.0);
    world.add(std::make_shared<Sphere>(Point3(4, 1, 0), 1.0, material3));
*/
    return world;
}

Color ray_color(const Ray& r, const Hittable& world, int depth) {
    hit_record rec;

    // If we've exceeded the ray bounce limit, no more light is gathered.
    if (depth <= 0)
        return Color(0,0,0);

    if (world.hit(r, 0.001, INF, rec)) {
        Ray scattered;
        Color attenuation;
        if (rec.mat_ptr->scatter(r, rec, attenuation, scattered))
            return attenuation * ray_color(scattered, world, depth-1);
        return Color(0,0,0);
    }
    Vec3 unit_direction = unit_vector(r.direction());
    auto t = 0.5*(unit_direction.y() + 1.0);
    return (1.0-t) * Color(1.0, 1.0,  1.0) + t * Color(0.5, 0.7, 1.0);
}

int main() {
    //19:50
    // Image

    const double aspect_ratio = 4.0/3.0;

    const int image_width = 500;
    const int image_height = static_cast<int>(image_width / aspect_ratio);
    const int samples_per_pixel = 20;
    const int max_depth = 50;
/*
    const double aspect_ratio = 16.0 /9.0;
    const int image_width = 256;
    const int image_height = static_cast<int>(image_width / aspect_ratio);
    const int samples_per_pixel = 100;
    const int max_depth = 50;
    // World
    HittableList world = random_scene();
    //Camera
    Vec3 vup(0,1,0);
    auto dist_to_focus = 10.0;
    Color background(0,0,0);
    Point3 lookfrom = Point3(12, 2, 1);
    Point3 lookat = Point3(0, 0, 0);
    auto aperture = 0.1;
    double vfov = 20.0;
*/

    HittableList world = cornell_box();
    //hitable *world = cornell_balls();
    //hitable *world = cornell_smoke();
    //hitable *world = cornell_final();
    //hitable *world = final();

    Vec3 lookfrom(278, 278, -800);
    //vec3 lookfrom(478, 278, -600);
    Vec3 lookat(278,278,0);
    //vec3 lookfrom(0, 0, 6);
    //vec3 lookat(0,0,0);
    float dist_to_focus = 10.0;
    float aperture = 0.0;
    float vfov = 40.0;

    Vec3 vup(0,1,0);
    Camera cam(lookfrom, lookat, vup, vfov, aspect_ratio, aperture, dist_to_focus);
    // Render
    if (MULTITHREAD){
        Color* pixelBuffer = (Color*) malloc(image_height * image_width * sizeof(Color));
        auto start = std::chrono::high_resolution_clock::now();
        omp_set_num_threads(16);
        #pragma omp parallel for
        for (int j = image_height-1; j >= 0; --j) {
            //std::cout << "\rScanned line: " << j << ' ' << std::endl;

            for (int i = 0; i < image_width; ++i) {
                Color pixel_color(0, 0, 0);
                for (int s = 0; s < samples_per_pixel; ++s) {
                    double u = (i + random_double()) / (image_width-1);
                    double v = (j + random_double()) / (image_height-1);
                    Ray r = cam.get_ray(u, v);
                    //#pragma omp parallel
                    pixel_color += ray_color(r, world, max_depth);
                }
                pixelBuffer[j*image_width+i] = pixel_color;
            }
        }
        auto stop = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
        std::cout << "\nRUN TIME: " <<duration.count()/1000000.0 << std::endl;

        std::ofstream output;
        output.open("output.ppm");
        output << "P3\n" << image_width << " " << image_height << "\n255\n";
        for (int j = image_height-1; j >= 0; --j) {
            for (int i = 0; i < image_width; ++i) {
                write_color(output, pixelBuffer[j*image_width+i], samples_per_pixel);
            }
        }
        output.close();
        free(pixelBuffer);
    }
    else{

        std::ofstream output;
        output.open("output.ppm");
        output << "P3\n" << image_width << " " << image_height << "\n255\n";
         for (int j = image_height-1; j >= 0; --j) {
             std::cerr << "\rScanlines remaining: " << j << ' ' << std::flush;
             for (int i = 0; i < image_width; ++i) {
                 Color pixel_color(0, 0, 0);
                 for (int s = 0; s < samples_per_pixel; ++s) {
                     double u = (i + random_double()) / (image_width-1);
                     double v = (j + random_double()) / (image_height-1);
                     Ray r = cam.get_ray(u, v);
                     pixel_color += ray_color(r, world, max_depth);
                 }
                 write_color(output, pixel_color, samples_per_pixel);
             }
         }

        output.close();

        std::cerr << "\nDone.\n";
    }

}