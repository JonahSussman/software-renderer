#include "sw_math.hpp"

#include <cmath>

#include <iostream>
#include <tuple>

using namespace std;
using namespace math;

int main() {
  /*
  math::Matrix<3, 4> m_0 = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
  math::Matrix<4, 3> m_1 = {11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0};
  auto m_2 = m_0 * m_1;

  double x = math::Matrix<1, 1>({ 1 });

  for (auto x : m_2) std::cout << x << " ";
  std::cout << std::endl;

  std::cout << "[";
  for (int row = 0; row < 3; row++) {
    for (int col = 0; col < 4; col++) {
      std::cout << " " << m_0(row, col);
      if (col != 4-1) std::cout << ",";
    } if (row != 3-1) std::cout << ";";
  }
  std::cout << " ]" << std::endl;

  math::IdentityMatrix<4> I4;
  math::Matrix<4, 4> i4 = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
  math::Matrix<4, 4> Ii44 = I4 + i4;

  math::Vec3 X = 3.0 * math::Vec3({1, 0, 0});

  std::cout << "[";
  for (int row = 0; row < 4; row++) {
    for (int col = 0; col < 4; col++) {
      std::cout << " " << I4(row, col);
      if (col != 4-1) std::cout << ",";
    } if (row != 3-1) std::cout << ";";
  }

  std::cout << " ]" << std::endl;

  for (auto x : Ii44.row(3)) {
    std::cout << x << " ";
  } std::cout << std::endl;

  math::Vec4 v0 = { 1, 2, 3, 4 };
  math::Vec4 v1 = { 4, 3, 2, 1 };

  std::cout << math::dot(v0, v1) << std::endl;

///////////////////////////////////////////////
  // THE IDENTITY MATRIX
  math::Matrix<4,4> thing = {
    1, 0, 0, 0,
    0, 1, 0, 0,
    0, 0, 1, 0,
    0, 0, 0, 1
  };

  // A COLUMN VECTOR (aka math::Matrix<4, 1>)
  math::Vec<4> thing_v = { 
    1.0, 
    3.0, 
    6.0, 
    -13.0 
  };

  math::Vec<4> a = thing_v;

  for (int i = 0; i < 100; i++) a = thing * a;

  std::cout << a[0] << " " << a[1] << " " << a[2] << " " << a[3] << std::endl;
  std::cout << thing_v[0] << " " << thing_v[1] << " " << thing_v[2] << " " << thing_v[3] << std::endl;
  */
  ///////////////////////////////
  // ROTOR TESTS               //
  ///////////////////////////////

  Vec3 pt = { 1, 0, 0 };

  BiVec3 plane = { 1, 0, 0 };
  plane.normalize();
  Matrix<4, 4> rm = {
    1, 0, 0, 0,
    0, 1, 0, 0,
    0, 0, 1, 0,
    0, 0, 0, 1
  };
  Rotor3 rtr(plane, plane.length());
  /* Works? 
  for (int i = 0; i < 1000; i++) {
    rm = rtr.to_matrix4() * rm;
    Vec4 p = rm * Vec4({ pt[0], pt[1], pt[2], 1 });
    cout << p[0] << " " << p[1] << " " << p[2] << endl;
  }

  Rotor3 rot({ 1, 0, 0, 0 });
  for (int i = 0; i < 100; i++) {
    rot = rtr.rotate(rot);
    Vec3 p = rot.to_matrix3() * pt;
    cout << p[0] << " " << p[1] << " " << p[2] << endl;
  }
  */

  Vec3 v = {
    1, 2, 3
  };

  Vec4 r_0 = { v[0], v[1], v[2], 0 };
  for (int i = 0; i < 100; i++) {
    rm = rtr.to_matrix4() * rm;
  }
  r_0 = rm * r_0;

  auto r_1 = v;
  for (int i = 0; i < 100; i++) r_1 = rtr.rotate(r_1);

  Rotor3 R;
  for (int i = 0; i < 100; i++) R *= rtr;
  auto r_2 = R.rotate(v);

  cout << r_0[0] << " " << r_0[1] << " " << r_0[2] << endl;
  cout << r_1[0] << " " << r_1[1] << " " << r_1[2] << endl;
  cout << r_2[0] << " " << r_2[1] << " " << r_2[2] << endl;
}