#pragma once

#include <cmath>

#include <array>
#include <type_traits>

namespace math {

double lerp(double f_0, double f_1, double x_0, double x_1, double x) {
  return f_0 + ((x - x_0) / (x_1 - x_0)) * (f_1 - f_0);
}

double lerp(double f_0, double f_1, double t) {
  return f_0 + t * (f_1 - f_0);
}

/**
 * Matrix class with N rows and M columns. Stores type T. Subclasses std::array.
 * 
 * The data is stored as one contiguous block of memory in an std::array in row-
 * major order (the 1st row is stored, followed by the 2nd, followed by the 3rd,
 * etc...)
 */
template <std::size_t N, std::size_t M>
struct Matrix : std::array<double, N*M> {
  static_assert(
    N > 0 && M > 0,
    "A math::Matrix must have a positive number of rows and columns"
  );

  operator double() const {
    static_assert(N == 1 && M == 1);
    return (*this)[0];
  }

  double& operator()(int loc)                { return (*this)[loc]; }
  double  operator()(int loc) const          { return (*this)[loc]; }
  double& operator()(int row, int col)       { return (*this)[M*row + col]; }
  double  operator()(int row, int col) const { return (*this)[M*row + col]; }

  constexpr auto row_size() { return N; }
  constexpr auto col_size() { return M; }

  constexpr auto is_square() { return N == M; }

  Matrix<1, M> row(int row) { 
    Matrix<1, M> result{};
    for (int col = 0; col < M; col++) result(col) = (*this)(row, col);
    return result;
  }

  Matrix<N, 1> col(int col) {
    Matrix<N, 1> result{};
    for (int row = 0; row < N; row++) result(row) = (*this)(row, col);
    return result;
  }

  Matrix<M, N> transpose() {
    Matrix<M, N> result{};
    for (int row = 0; row < N; row++)
      for (int col = 0; col < M; col++)
        result(col, row) = (*this)(row, col);

    return result;
  }

  double length_squared() const {
    static_assert(
      N == 1 || M == 1, 
      "Matrix::length_squared() is only valid for 1xM or Nx1 matrices."
    );

    double result = 0.0;
    for (int i = 0; i < N*M; i++) result += (*this)[i] * (*this)[i];
    return result;
  }

  double length() const {
    static_assert(
      N == 1 || M == 1, 
      "Matrix::length() is only valid for 1xM or Nx1 matrices."
    );

    return sqrt(length_squared());
  }

  void normalize() {
    static_assert(
      N == 1 || M == 1, 
      "Matrix::normalize() is only valid for 1xM or Nx1 matrices."
    );

    double l = length();
    for (int i = 0; i < N*M; i++) (*this)[i] /= l;
  }

  Matrix<N, M> normal() {
    static_assert(
      N == 1 || M == 1, 
      "Matrix::normalize() is only valid for 1xM or Nx1 matrices."
    );

    Matrix<N, M> result = (*this);
    result.normalize();
    return result;
  }

  Matrix<N, M> operator+(const Matrix<N, M>& other) const {
    Matrix<N, M> result = *this;
    for (int i = 0; i < N*M; i++) result[i] += other[i];
    return result;
  }

  Matrix<N, M> operator+=(const Matrix<N, M>& other) {
    for (int i = 0; i < N*M; i++) (*this)[i] += other[i];
    return (*this);
  }

  Matrix<N, M> operator-(const Matrix<N, M>& other) const {
    Matrix<N, M> result = *this;
    for (int i = 0; i < N*M; i++) result[i] -= other[i];
    return result;
  }

  Matrix<N, M> operator*(const double& scalar) {
    Matrix<N, M> result = *this;
    for (int i = 0; i < N*M; i++) result[i] *= scalar;
    return result;
  }

  friend Matrix<N, M> operator*(const double& scalar, const Matrix<N, M>& m) {
    Matrix<N, M> result = m;
    for (int i = 0; i < N*M; i++) result[i] *= scalar;
    return result;
  }

  template <std::size_t K>
  auto operator*(const Matrix<M, K>& other) {
    Matrix<N, K> result{};

    for (int row = 0; row < N; row++) {
      for (int col = 0; col < K; col++) {
        for (int term = 0; term < M; term++) {
          result(row, col) += (*this)(row, term) * other(term, col);
        }
      }
    }

    if constexpr (N == 1 && K == 1) {
      return (double)result;
    } else {
      return result;
    } 
  }
};

/**
 * Creates an NxN identity matrix. Subclasses Matrix.
 */
template <std::size_t N>
struct IdentityMatrix : Matrix<N, N> {
  IdentityMatrix() {
    for (int i = 0; i < N; i++) 
      for (int j = 0; j < N; j++)
        (*this)(i, j) = i == j ? 1.0 : 0.0;
  }
};

template<std::size_t N>
using Vec    = Matrix<N, 1>;

using Vec3   = Matrix<3, 1>;
using BiVec3 = Matrix<3, 1>;
using Vec4   = Matrix<4, 1>;

/**
 * Computes the dot product of two std::arrays. (!!) Because matrix subclasses
 * std::array, and a 'vector' is a math::matrix, this is valid and removes all
 * the weird edge cases. Though if applied wrong, it can have some funny results
 * such as summing up the element-wise multiplication of two matrices.
 */
template <std::size_t N>
double dot(std::array<double, N> left, std::array<double, N> right) {
  double result = 0.0;
  for (int i = 0; i < N; i++) result += left[i] * right[i];
  return result;
}

/**
 * Computes the wedge product (from geometric algebra) of two vectors.
 */
BiVec3 wedge(const Vec3& u, const Vec3& v) {
  BiVec3 result = {
    u[0]*v[1] - u[1]*v[0], // XY
    u[0]*v[2] - u[2]*v[0], // XZ
    u[1]*v[2] - u[2]*v[1]  // YZ
  };
  return result;
}

/**
 * This math library uses rotors from Geometric Algebra to handle rotation (as
 * opposed to quaternions).
 */
struct Rotor3 {
  double a; 
  BiVec3 b;

  Rotor3() {
    a = 1.0;
    b = { 0.0, 0.0, 0.0 };
   }
  Rotor3(double a, double xy, double xz, double yz) :
    a(a), b({ xy, xz, yz }) { }
  
  Rotor3(Vec4 v) : a(v[0]), b({ v[1], v[2], v[3] }) { }
  Rotor3(double a, BiVec3 b) : a(a), b(b) { }

  Rotor3(const Vec3& v_from, const Vec3& v_to) {
    a = 1 + dot(v_to, v_from);
    b = wedge(v_to, v_from);
    normalize();
  }

  // Plane must be normalized
  Rotor3(const BiVec3& bv_plane, double angle_radian) {
    double sin_a = sin(angle_radian / 2.0);
    a = cos(angle_radian / 2.0);
    b[0] = -sin_a * bv_plane[0];
    b[1] = -sin_a * bv_plane[1];
    b[2] = -sin_a * bv_plane[2];
    // normalize();
  }

  Rotor3 operator*(const Rotor3& q) const {
    const Rotor3& p = *this;
    Rotor3 r;

    // r.a = p.a * q.a - p.b[0] * p.b[0] - p.b[1] * p.b[1] - p.b[2] * p.b[2];
    // r.b = {
    //   p.b[0] * q.a + p.a * q.b[0] + p.b[2] * q.b[1] - p.b[1] * q.b[2],
    //   p.b[1] * q.a + p.a * q.b[1] - p.b[2] * q.b[0] + p.b[0] * q.b[2],
    //   p.b[2] * q.a + p.a * q.b[2] + p.b[1] * q.b[0] - p.b[0] * q.b[1]
    // };

    r.a =    p.a    * q.a - p.b[0] * q.b[0] - p.b[1] * q.b[1] - p.b[2] * q.b[2];
    r.b[0] = p.b[0] * q.a + p.a    * q.b[0] + p.b[2] * q.b[1] - p.b[1] * q.b[2];
    r.b[1] = p.b[1] * q.a + p.a    * q.b[1] - p.b[2] * q.b[0] + p.b[0] * q.b[2];
    r.b[2] = p.b[2] * q.a + p.a    * q.b[2] + p.b[1] * q.b[0] - p.b[0] * q.b[1];

    return r;
  }

  Vec3 rotate(const Vec3& x) const {
    const Rotor3& p = *this;
    Vec3 q = {
      p.a * x[0] + x[1] * p.b[0] + x[2] * p.b[1],
      p.a * x[1] - x[0] * p.b[0] + x[2] * p.b[2],
      p.a * x[2] - x[0] * p.b[1] - x[1] * p.b[2]
    };

    double q012 = x[0] * p.b[2] - x[1] * p.b[1] + x[2] * p.b[0];

    Vec3 r = {
      p.a * q[0] + q[1] * p.b[0] + q[2] * p.b[1] + q012 * p.b[2],
      p.a * q[1] - q[0] * p.b[0] - q012 * p.b[1] + q[2] * p.b[2],
      p.a * q[2] + q012 * p.b[0] - q[0] * p.b[1] - q[1] * p.b[2]
    };

    return r;
  }

  Rotor3& operator*=(const Rotor3& r) {
    (*this) = (*this) * r;
    return *this;
  }

  Rotor3 rotate(const Rotor3& r) const {
    return (*this) * r * (*this).reverse();
  }

  Rotor3 reverse() const {
    return Rotor3(a, -1.0 * b[0], -1.0 * b[1], -1.0 * b[2]);
  }

  double length_squared() const {
    // return a*a + b.length_squared();
    return a*a + b[0]*b[0] + b[1]*b[1] + b[2]*b[2];
  }

  double length() const {
    return sqrt(length_squared());
  }

  void normalize() {
    double l = length();
    a /= l, b[0] /= l, b[1] /= l, b[2] /= l;
  }

  // Normalized rotor
  Rotor3 normal() const {
    Rotor3 r = *this;
    r.normalize();
    return r;
  }

  Matrix<3, 3> to_matrix3() const {
    Vec3 i = rotate({ 1.0, 0.0, 0.0 });
    Vec3 j = rotate({ 0.0, 1.0, 0.0 });
    Vec3 k = rotate({ 0.0, 0.0, 1.0 });

    return {
      i[0], j[0], k[0],
      i[1], j[1], k[1],
      i[2], j[2], k[2]
    };
  }

  // Homogenous coordinates
  Matrix<4, 4> to_matrix4() const {
    Vec3 i = rotate({ 1.0, 0.0, 0.0 });
    Vec3 j = rotate({ 0.0, 1.0, 0.0 });
    Vec3 k = rotate({ 0.0, 0.0, 1.0 });

    return {
      i[0], j[0], k[0], 0.0,
      i[1], j[1], k[1], 0.0,
      i[2], j[2], k[2], 0.0, 
      0.00, 0.00, 0.00, 1.0
    };
  }
};

Rotor3 geo(const Vec3& a, const Vec3& b) {
  return Rotor3(dot(a, b), wedge(a, b));
}

};