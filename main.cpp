#include "render.h"
#include <vector>
#include <cstdio>

int32_t main() {
  std::mt19937_64 mt_rand(std::chrono::system_clock::now().time_since_epoch().count());
  std::uniform_real_distribution<long double> urd(0, 1);

  std::vector<Sphere> spheres(9);

  spheres[0] = Sphere(1e5, Vector(1e5 + 1, 40.8, 81.6), Vector(), Vector(.75, .25, .25), DIFFuse);//Left
  spheres[1] = Sphere(1e5, Vector(-1e5 + 99, 40.8, 81.6), Vector(), Vector(.25, .25, .75), DIFFuse);//Right
  spheres[2] = Sphere(1e5, Vector(50, 40.8, 1e5), Vector(), Vector(.75, .75, .75), DIFFuse);//Back
  spheres[3] = Sphere(1e5, Vector(50, 40.8, -1e5 + 170), Vector(), Vector(), DIFFuse);//Front
  spheres[4] = Sphere(1e5, Vector(50, 1e5, 81.6), Vector(), Vector(.75, .75, .75), DIFFuse);//Bottom
  spheres[5] = Sphere(1e5, Vector(50, -1e5 + 81.6, 81.6), Vector(), Vector(.75, .75, .75), DIFFuse);//Top
  spheres[6] = Sphere(16.5, Vector(27, 16.5, 47), Vector(), Vector(0.99, 0.99, 0.99), SPECular);//Mirror
  spheres[7] = Sphere(16.5, Vector(73, 16.5, 78), Vector(), Vector(0.99, 0.99, 0.99), REFRactive);//Glass
  spheres[8] = Sphere(600, Vector(50, 681.6 - .27, 81.6), Vector(12, 12, 12), Vector(), DIFFuse);//Light

  size_t w = 1024, h = 768, samps = 5000 / 4;
  Vector cam_o(50, 52, 295.6);
  Vector cam_d(Vector(0, -0.042612, -1).norm());
  Ray cam(cam_o, cam_d);
  Vector cx = Vector(w * 0.5135L / h, 0, 0);
  Vector cy = (cx % cam.dir_).norm() * 0.5135L;
  Vector r{};

  std::vector<Vector> c(w * h);

#pragma omp parallel for schedule(dynamic, 1) private(r)
  {
    for (size_t y(0); y < h; y++) {
      printf("\rRendering (%lu spp) %5.2Lf%%...", samps * 4, 100.0L * y / (h - 1));
      for (size_t x = 0; x < w; x++)
        for (size_t sy = 0, i = (h - y - 1) * w + x; sy < 2; sy++)
          for (size_t sx = 0; sx < 2; sx++, r = Vector()) {
            for (size_t s = 0; s < samps; s++) {
              long double r1 = 2 * urd(mt_rand), dx = r1 < 1 ? sqrtl(r1) - 1 : 1 - sqrtl(2 - r1);
              long double r2 = 2 * urd(mt_rand), dy = r2 < 1 ? sqrtl(r2) - 1 : 1 - sqrtl(2 - r2);
              Vector d = cam.dir_ + \
                  cx * (((sx + 0.5L + dx) / 2 + x) / w - 0.5L) + \
                  cy * (((sy + 0.5L + dy) / 2 + y) / h - 0.5L);
              r = r + trace(spheres, Ray(cam.o_ + d * 140, d.norm()), 0) * (1.0L / samps);
            }
            c[i] = c[i] + Vector(clamp(r.x_), clamp(r.y_), clamp(r.z_)) * 0.25L;
          }
    }
  }

  printf("Writing Image...\n");
  FILE *f = fopen("image.ppm", "w");
  fprintf(f, "P3\n%lu %lu\n%d\n", w, h, 255);
  for (int i = 0; i < w * h; i++)
    fprintf(f, "%d %d %d ", toInt(c[i].x_), toInt(c[i].y_), toInt(c[i].z_));
  fclose(f);
  return 0;
}