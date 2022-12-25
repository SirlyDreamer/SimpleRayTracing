#include <cmath>

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
 private:
  long double x_, y_, z_;
};

class Ray {
 public:
  Ray() = default;
  Ray(Vector &o, Vector &dir) : o_(o), dir_(dir) {}
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
 private:
  static constexpr long double eps = 1e-5;
  long double r_;
  Vector o_, dir_, color_;
  refl_t refl_;
};