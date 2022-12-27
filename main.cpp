#include "render.h"
#include <vector>
#include <cstdio>

void initSphere(std::vector<Sphere> &sp) {
  Vector Cen(50, 40.8, -860);
  sp.emplace_back(1600,
                  Vector(1, 0, 2) * 3000,
                  Vector(1, .9, .8) * 1.2e1 * 1.56 * 2,
                  Vector(),
                  DIFFuse);
  sp.emplace_back(1560,
                  Vector(1, 0, 2) * 3500,
                  Vector(1, .5, .05) * 4.8e1 * 1.56 * 2,
                  Vector(),
                  DIFFuse);
  sp.emplace_back(10000,
                  Cen + Vector(0, 0, -200),
                  Vector(0.0627, 0.188, 0.569) * 6e-2 * 8,
                  Vector(.7, .7, 1) * .25,
                  DIFFuse);
  sp.emplace_back(10000,
                  Cen + Vector(0, 0, -200),
                  Vector(0.00063842, 0.02001478, 0.28923243) * 6e-2 * 8,
                  Vector(.7, .7, 1) * .25,
                  DIFFuse);
  sp.emplace_back(100000,
                  Vector(50, -100000, 0),
                  Vector(),
                  Vector(.3, .3, .3),
                  DIFFuse);
  sp.emplace_back(110000,
                  Vector(50, -110048.5, 0),
                  Vector(.9, .5, .05) * 4,
                  Vector(),
                  DIFFuse);
  sp.emplace_back(4e4,
                  Vector(50, -4e4 - 30, -3000),
                  Vector(),
                  Vector(.2, .2, .2),
                  DIFFuse);
  sp.emplace_back(3.99e4,
                  Vector(50, -3.99e4 + 20.045, -3000),
                  Vector(),
                  Vector(.7, .7, .7),
                  DIFFuse);
  sp.emplace_back(26.5,
                  Vector(22, 26.5, 42),
                  Vector(),
                  Vector(1, 1, 1) * .596,
                  SPECular);
  sp.emplace_back(13,
                  Vector(75, 13, 82),
                  Vector(),
                  Vector(.96, .96, .96) * .96,
                  REFRactive);
  sp.emplace_back(22,
                  Vector(87, 22, 24),
                  Vector(),
                  Vector(.6, .6, .6) * .696,
                  REFRactive);
}

int32_t main() {
  std::mt19937_64 mt_rand(std::chrono::system_clock::now().time_since_epoch().count());
  std::uniform_real_distribution<long double> urd(0, 1);

  std::vector<Sphere> spheres;
  initSphere(spheres);

  size_t w = 1920, h = 1080, samps = 25000 / 4;
  Vector cam_o(50, 52, 295.6);
  Vector cam_d(Vector(0, -0.042612, -1).norm());
  Ray cam(cam_o, cam_d);
  Vector cx = Vector(w * 0.5135L / h, 0, 0);
  Vector cy = (cx % cam.dir_).norm() * 0.5135L;
  Vector r{};

  std::vector<Vector> c(w * h);

#pragma omp parallel for schedule(dynamic, 1) private(r)
  for (size_t y = 0; y < h; y++) {
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

  printf("Writing Image...\n");
  FILE *f = fopen("image.ppm", "w");
  fprintf(f, "P3\n%lu %lu\n%d\n", w, h, 255);
  for (int i = 0; i < w * h; i++)
    fprintf(f, "%d %d %d ", toInt(c[i].x_), toInt(c[i].y_), toInt(c[i].z_));
  fclose(f);
  return 0;
}