#pragma once

#include <cmath>

#include <array>
#include <iostream>
#include <iomanip>

/**
 * Condensed linear algebra library for 4x4 matrices.
 */
namespace fmth {

/**
 * 4x1 column vector (technically it could also be a row vector).
 */
struct Vec {
  double x, y, z, w;

  double length_squared() { return x*x + y*y + z*z + w*w; }
  double length()         { return sqrt(length_squared()); }

  Vec normalized() { auto l = length(); return { x/l, y/l, z/l, w/l }; }
  
  Vec  operator-()              const { return {    -x,    -y,    -z,    -w }; }
  Vec  operator+ (const Vec& v) const { return { x+v.x, y+v.y, z+v.z, w+v.w }; }
  Vec  operator- (const Vec& v) const { return { x-v.x, y-v.y, z-v.z, w-v.w }; }
  Vec  operator* (double s)     const { return {   s*x,   s*y,   s*z,   s*w }; }
  Vec& operator+=(const Vec& v)       { *this = (*this) + v; return *this; }
  Vec& operator-=(const Vec& v)       { *this = (*this) - v; return *this; }
  Vec& operator*=(double s)           { *this = (*this) * s; return *this; }
  friend Vec operator*(double s, const Vec& v) { return v * s; }

  // Shifts all the indicies left by x, replacing them with zeros.
  Vec  operator<<(int i) const {
    std::array<double, 8> a = { x, y, z, w, 0, 0, 0, 0 };
    if (i > 4) i = 4; if (i < 0) i = 0;
    return { a[i], a[i+1], a[i+2], a[i+3] };
  }

  // Shifts all the indicies right by x, replacing them with zeros.
  Vec  operator>>(int i) const {
    std::array<double, 8> a = { 0, 0, 0, 0, x, y, z, w };
    if (i > 4) i = 4; if (i < 0) i = 0;
    return { a[4-i], a[5-i], a[6-i], a[7-i] };
  }
};

/**
 * Returns dot product of two 4x1 vectors
 */
double dot(const Vec& a, const Vec& b) {
  return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w;
}

/**
 * MEGA HYPER ABUSE OF MATHEMATICS ALERT! 
 * Returns the "wedge" product of two 3x1 vectors. However the vector inputs are
 * stored as 4x1 Vecs and the output is ALSO a 4x1 Vec. Technically, what this 
 * actually is is the cross product. For the true wedge product the input should
 * be a Vec3 and the output should be a BiVec3, but this way of doing it just
 * makes implementation easier.
 */
Vec wedge(const Vec& u, const Vec& v) {
  return {
    u.x*v.y - u.y*v.x, // xy plane
    u.x*v.z - u.z*v.x, // xz plane
    u.y*v.z - u.z*v.y, // yz plane
    (u.w + v.w) / 2.0  // Homogenous coordinates wizardry
  };
}

/**
 * 4x4 Matrix class.
 * NOTE: The matrix is stored as it's transposition, i.e. if the column vectors
 * are u, v, w and x then it is stored as { u.x, u.y, u.z, u.w, v.x, ..., x.x }.
 * This is so the column vectors are fast and easy to remove, splice, etc...
 */
struct Mat {
  std::array<double, 16> data;

  // Default ctors.
  Mat() : data({}) { }
  Mat(const Mat& m) : data(m.data) { }
  Mat(std::array<double, 16> list) : data(list) { }
  Mat& operator=(const Mat& m) { if (this != &m) data = m.data; return *this; }
  
  // Construct a matrix out of columns
  Mat(Vec c0, Vec c1, Vec c2, Vec c3) {
    data = {
      c0.x, c0.y, c0.z, c0.w,
      c1.x, c1.y, c1.z, c1.w,
      c2.x, c2.y, c2.z, c2.w,
      c3.x, c3.y, c3.z, c3.w
    };
  }

  inline double& at(int row, int col)       { return data[4*col + row]; }
  inline double  at(int row, int col) const { return data[4*col + row]; }

  // Abuse the fact structs are stored sequentially and copy the memory
  Vec col(int c) const { 
    Vec result;
    memcpy(&result, data.data() + 4*c, sizeof(Vec));
    return result;
  }

  Vec row(int r) const { return { data[r], data[r+4], data[r+8], data[r+12] }; }

  std::array<Vec, 4> cols() const { return { col(0), col(1), col(2), col(3) }; }
  std::array<Vec, 4> rows() const { return { row(0), row(1), row(2), row(3) }; }

  inline void set_col(int c, const Vec v) { 
    memcpy(data.data() + 4*c, &v, sizeof(Vec));
  }

  Mat transpose() const {
    Mat m = *this;

    for (int row = 0; row < 4; ++row) {
      for (int col = 0; col < 4; ++col) {
        m.at(row, col) = (*this).at(col, row);
      }
    }

    return m;
  }

  Mat operator*(const Mat& m) const {
    Mat result;

    for (int row = 0; row < 4; ++row) {
      for (int col = 0; col < 4; ++col) {
        for (int term = 0; term < 4; ++term) {
          result.at(row, col) += (*this).at(row, term) * m.at(term, col);
        }
      }
    }

    return result;
  }

  Vec operator*(const Vec& v) const {
    auto r = rows();
    return { dot(v, r[0]), dot(v, r[1]), dot(v, r[2]), dot(v, r[3]) };
  }
};

const Mat IdentityMat = Mat({
  1, 0, 0, 0,
  0, 1, 0, 0,
  0, 0, 1, 0,
  0, 0, 0, 1
});

struct Rot {
  Vec r;

  Rot() : r({ 1, 0, 0, 0 }) { }
  Rot(const Vec& v) : r(v) { }
  Rot(double a, double xy, double xz, double yz) : r({ a, xy, xz, yz }) { }
  
  Rot(const Vec& from, const Vec& to) {
    Vec a = { 1 + dot(to, from), 0, 0, 0 };
    Vec b = wedge(from, to) >> 1;
    r = (a + b).normalized();
  }

  // Plane must be normalized
  Rot(const Vec& plane, double angle) {
    double s = sin(angle / 2.0);
    r = { cos(angle / 2.0), -s*plane.x, -s*plane.y, -s*plane.z };
  }

  Rot reverse() const { return Rot(r.x, -r.y, -r.z, -r.w ); }
  void normalize() { r = r.normalized(); }

  Vec rotate(const Vec& v) const {
    const Rot& p = *this;
    Vec q = {
      p.r.x * v.x + v.y * p.r.y + v.z * p.r.z,
      p.r.x * v.y - v.x * p.r.y + v.z * p.r.w,
      p.r.x * v.z - v.x * p.r.z - v.y * p.r.w,
      v.w
    };

    double xyz = v.x * p.r.w - v.y * p.r.z + v.z * p.r.y;

    Vec r = {
      p.r.x * q.x + q.y * p.r.y + q.z * p.r.z + xyz * p.r.w,
      p.r.x * q.y - q.x * p.r.y - xyz * p.r.z + q.z * p.r.w,
      p.r.x * q.z + xyz * p.r.y - q.x * p.r.z - q.y * p.r.w,
      q.w
    };

    return r;
  }

  Mat matrix() const {
    Vec i = rotate({ 1, 0, 0, 0 });
    Vec j = rotate({ 0, 1, 0, 0 });
    Vec k = rotate({ 0, 0, 1, 0 });
    Vec l =        { 0, 0, 0, 1 };

    return Mat(i, j, k, l);
  }
};

std::ostream& operator<<(std::ostream& o, const Vec& v) {
  o << "[ " << std::setw(8) << v.x 
    <<  " " << std::setw(8) << v.y 
    <<  " " << std::setw(8) << v.z 
    <<  " " << std::setw(8) << v.w << " ]";
}

std::ostream& operator<<(std::ostream& o, const Mat& m) {
  for (int row = 0; row < 4; ++row) {
    o << "[ ";
    for (int col = 0; col < 4; ++col) {
      o << std::setw(8) << m.at(row, col) << " ";
    } o << (row == 3 ? "]" : "]\n");
  }
}

}