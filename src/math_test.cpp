#include "sw_math.hpp"

#include <iostream>

int main() {
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
}