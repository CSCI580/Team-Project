#include "util.h"
#include "Color.h"
#include "Hittable.h"
#include "HittableList.h"
#include "Sphere.h"
#include "Camera.h"

#include <memory>
#include <iostream>
#include <fstream>
#include <random>

bool valid_bubble(Vec3 bubble_offset, double bubble_size) {
    if (bubble_size < 0 || bubble_size > 1) {
        return false;
    }
    for (int i = 0; i < 3; ++i) {
        if (bubble_offset[i] < -1 + bubble_size || bubble_offset[i] > 1 - bubble_size) {
            return false;
        }
    }
    return true;
}

void generate_bubble(HittableList* world, Vec3 center) {
    const int BUBBLE_NUM = 300;
    std::default_random_engine generator;
    std::normal_distribution<double> pos_dist(0, 0.2);
    std::normal_distribution<double> size_dist(0.02,0.003);

    auto bubble_material = std::make_shared<dielectric>(1/1.31);
    for (int i = 0; i < BUBBLE_NUM; ++i) {
        Vec3 bubble_offset;
        double bubble_size;
        do {
            bubble_offset = Vec3(pos_dist(generator), pos_dist(generator), pos_dist(generator));
            bubble_size = size_dist(generator);
        } while (!valid_bubble(bubble_offset, bubble_size));
        std::cout << bubble_offset+center << "   " << bubble_size << "\n";
        world->add(std::make_shared<Sphere>(bubble_offset+center, bubble_size, bubble_material));
    }
}

HittableList dielectric_scene() {
    HittableList world;

    auto ground_material = std::make_shared<lambertian>(Color(0.5, 0.5, 0.5));
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
                    // metal
                    auto albedo = Color::random(0.5, 1);
                    auto fuzz = random_double(0, 0.5);
                    sphere_material = std::make_shared<metal>(albedo, fuzz);
                    world.add(std::make_shared<Sphere>(center, 0.2, sphere_material));
                    // glass
//                    sphere_material = std::make_shared<dielectric>(1.5);
//                    world.add(std::make_shared<Sphere>(center, 0.2, sphere_material));
                }
            }
        }
    }
//    auto material0 = std::make_shared<dielectric>(1);
//    world.add(std::make_shared<Sphere>(Point3(0, 1, 0), 1.01, material0));

    auto material1 = std::make_shared<dielectric>(1.31);
    world.add(std::make_shared<Sphere>(Point3(0, 1, 0), 1.0, material1));

//    auto bubble_material = std::make_shared<dielectric>(1/1.31);
//    world.add(std::make_shared<Sphere>(Point3(0.205287, 1-0.326045, -0.0365897), 0.5, bubble_material));


    generate_bubble(&world, Vec3(0, 1, 0));

//
//    auto material3 = std::make_shared<metal>(Color(0.7, 0.6, 0.5), 0.0);
//    world.add(std::make_shared<Sphere>(Point3(4, 1, 0), 1.0, material3));

    return world;
}

HittableList random_scene() {
    HittableList world;

    auto ground_material = std::make_shared<lambertian>(Color(0.5, 0.5, 0.5));
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

    auto material1 = std::make_shared<dielectric>(1.31);
    world.add(std::make_shared<Sphere>(Point3(0, 1, 0), 1.0, material1));

    auto material2 = std::make_shared<lambertian>(Color(0.4, 0.2, 0.1));
    world.add(std::make_shared<Sphere>(Point3(-4, 1, 0), 1.0, material2));

    auto material3 = std::make_shared<metal>(Color(0.7, 0.6, 0.5), 0.0);
    world.add(std::make_shared<Sphere>(Point3(4, 1, 0), 1.0, material3));

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
    const double aspect_ratio = 16.0 / 9.0;
    const int image_width = 512;
    const int image_height = static_cast<int>(image_width / aspect_ratio);
    const int samples_per_pixel = 30;
    const int max_depth = 50;

    // World
    HittableList world = dielectric_scene();

    //Camera
    Point3 lookfrom(13,2,3);
    Point3 lookat(0,0,0);
    Vec3 vup(0,1,0);
    auto dist_to_focus = 10.0;
    auto aperture = 0.1;

    Camera cam(lookfrom, lookat, vup, 20, aspect_ratio, aperture, dist_to_focus);
    // Render
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