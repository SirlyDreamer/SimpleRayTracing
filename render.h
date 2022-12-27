#include <cmath>
#include <chrono>
#include <random>
#include <vector>
#include <algorithm>

class Vector {
 public:
  Vector() = default;

  Vector(long double x, long double y, long double z) : x_(x), y_(y), z_(z) {}

  Vector operator+(const Vector &e) const {
    return Vector(x_ + e.x_, y_ + e.y_, z_ + e.z_);
  }

  Vector operator-(const Vector &e) const {
    return Vector(x_ - e.x_, y_ - e.y_, z_ - e.z_);
  }

  Vector operator*(long double lambda) const {
    return Vector(x_ * lambda, y_ * lambda, z_ * lambda);
  }

  Vector operator*(const Vector &e) const {
    return Vector(x_ * e.x_, y_ * e.y_, z_ * e.z_);
  }

  Vector &norm() {
    return *this = *this * (1 / sqrtl(x_ * x_ + y_ * y_ + z_ * z_));
  }

  long double dot(const Vector &e) const {
    return x_ * e.x_ + y_ * e.y_ + z_ * e.z_;
  }

  Vector operator%(const Vector &e) const {
    return Vector(y_ * e.z_ - z_ * e.y_, z_ * e.x_ - x_ * e.z_, x_ * e.y_ - y_ * e.x_);
  }

  long double x_, y_, z_;
};

class Ray {
 public:
  Ray() = default;

  Ray(Vector o, Vector dir) : o_(o), dir_(dir) {}

  Vector o_, dir_;
};

enum refl_t { DIFFuse, SPECular, REFRactive };

class Sphere {
 public:
  Sphere() = default;

  Sphere(long double r, Vector o, Vector dir, Vector color, refl_t refl) :
      r_(r), o_(o), dir_(dir), color_(color), refl_(refl) {}

  long double intersect(const Ray &ray) const {
    // d·d * t^2 + 2*(o-p)·d * t + (o-p)·(o-p)-R^2 = 0
    Vector o_p = o_ - ray.o_;
    long double b = o_p.dot(ray.dir_);
    long double delta = b * b - o_p.dot(o_p) + r_ * r_;
    if (delta < 0)
      return 0;
    delta = sqrtl(delta);
    long double t1 = b - delta, t2 = b + delta;
    if (t1 > eps)
      return t1;
    if (t2 > eps)
      return t2;
    return 0;
  }

  static constexpr long double eps = 1e-4;
  long double r_;
  Vector o_, dir_, color_;
  refl_t refl_;
};

inline long double clamp(long double x) {
  if (x < 0)
    return 0;
  if (x > 1)
    return 1;
  return x;
}

inline int32_t toInt(long double x) {
  return powl(clamp(x), 1 / 2.2) * 255 + 0.5;
}

inline bool intersect(const std::vector<Sphere> &s, const Ray &r, long double &t, size_t &id) {
  constexpr long double inf = 1e20, eps = 1e-5;
  t = inf;
  for (int i(s.size() - 1); i >= 0; --i) {
    long double d = s[i].intersect(r);
    if (d != 0 && d < t) {
      t = d;
      id = i;
    }
  }
  return t < inf;
}

Vector trace(const std::vector<Sphere> &s, const Ray &r, int depth) {
  std::mt19937_64 mt_rand(std::chrono::system_clock::now().time_since_epoch().count());
  std::uniform_real_distribution<long double> urd(0, 1);
  long double t;
  size_t id(0);

  if (!intersect(s, r, t, id)) {
    return Vector{};
  }

  const Sphere &obj = s[id];
  Vector x = r.o_ + r.dir_ * t;
  Vector n = (x - obj.o_).norm();
  Vector nl = n.dot(r.dir_) < 0 ? n : n * -1;
  Vector f = obj.color_;

  long double p = std::max({f.x_, f.y_, f.z_});

  if (depth > 4) {
    if (urd(mt_rand) < p)
      f = f * (1 / p);
    else
      return obj.dir_;
  }

  switch (obj.refl_) {
    case DIFFuse: {
      long double r1 = 2 * M_PI * urd(mt_rand);
      long double r2 = urd(mt_rand);
      long double r2s = sqrtl(r2);
      Vector w = nl;
      Vector u = ((std::abs(w.x_) > 0.1 ? Vector(0, 1, 0) : Vector(1, 0, 0)) % w).norm();
      Vector v = w % u;
      Vector d = (u * std::cos(r1) * r2s + v * std::sin(r1) * r2s + w * (sqrtl(1 - r2))).norm();
      Ray outRay(x, d);
      return obj.dir_ + f * trace(s, outRay, depth + 1);
    }
    case SPECular: {
      Ray reflRay(x, r.dir_ - n * 2 * n.dot(r.dir_));
      return obj.dir_ + f * trace(s, reflRay, depth + 1);
    }
    case REFRactive: {
      Ray reflRay(x, r.dir_ - n * 2 * n.dot(r.dir_));
      bool into = n.dot(nl) > 0;
      long double nc(1), nt(1.5);
      long double nnt = into ? nc / nt : nt / nc;
      long double ddn = r.dir_.dot(nl);
      long double cos2t = 1 - nnt * nnt * (1 - ddn * ddn);
      if (cos2t < 0)
        return obj.dir_ + f * trace(s, reflRay, depth + 1);
      Vector tdir = (r.dir_ * nnt - n * ((into ? 1 : -1) * (ddn * nnt + sqrtl(cos2t)))).norm();
      Ray outRay(x, tdir);
      long double a = nt - nc, b = nt + nc, R0 = a * a / (b * b);
      long double c = 1 - (into ? -ddn : tdir.dot(n));
      long double Re = R0 + (1 - R0) * c * c * c * c * c;
      long double Tr = 1 - Re;
      long double P = 0.25 + 0.5 * Re;
      long double RP = Re / p;
      long double TP = Tr / (1 - P);
      if (depth <= 1)
        return obj.dir_ + f * (trace(s, reflRay, depth + 1) * Re + trace(s, outRay, depth + 1) * Tr);
      if (urd(mt_rand) < P)
        return obj.dir_ + f * trace(s, reflRay, depth + 1) * RP;
      return obj.dir_ + f * trace(s, outRay, depth + 1) * TP;
    }
  }
}