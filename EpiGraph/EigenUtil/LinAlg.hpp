//
// Created by roryh on 19/09/2020.
//

#ifndef EPIGRAPH_EIGEN_UTIL_H
#define EPIGRAPH_EIGEN_UTIL_H

#include <EpiGraph/EigenUtil/StaticAsserts.hpp>

#include <Eigen/Dense>

#include <fstream>
#include <map>

namespace EpiGraph {

template <IsColVector Derived, typename F>
auto general_outer_product(const Derived &v, const Derived &u, F f)
    -> Eigen::Matrix<typename Derived::Scalar, Derived::RowsAtCompileTime,
                     Derived::RowsAtCompileTime> {

  using Scalar = typename Derived::Scalar;
  int rows = Derived::RowsAtCompileTime;

  Eigen::Matrix<Scalar, rows, rows> mat(v.rows(), u.rows());
  for (int j = 0; j < v.rows(); j++) {
    for (int i = 0; i < u.rows(); i++) {
      mat(i, j) = f(v[i], u[j]);
    }
  }
  return mat;
}

} // namespace EpiGraph
#endif // EPIGRAPH_EIGEN_UTIL_H
