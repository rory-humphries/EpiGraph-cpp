//
// Created by roryh on 22/06/2020.
//

#ifndef EPIGRAPH_SPATIAL_UTIL_H
#define EPIGRAPH_SPATIAL_UTIL_H

#include <Eigen/Dense>

namespace EpiGraph {
auto haversine(double theta) -> double;

auto long_lat_distance(double lon1, double lat1, double lon2, double lat2)
    -> double;

auto long_lat_distance_2(double lon1, double lat1, double lon2, double lat2)
    -> double;

template <typename Mat1, typename Mat2>
auto distance_matrix(Eigen::MatrixBase<Mat2> &pos_mat) -> Mat1 {
  static_assert(Mat2::ColsAtCompileTime == 2, "Expected a N x 2 matrix");

  Mat1 mat(pos_mat.rows(), pos_mat.rows());
  mat.setZero();

  for (int i = 0; i < pos_mat.rows(); i++) {
    for (int j = 0; j < pos_mat.rows(); j++) {

      double lon1 = pos_mat(i, 0), lon2 = pos_mat(j, 0), lat1 = pos_mat(i, 1),
             lat2 = pos_mat(j, 1);
      mat(i, j) = long_lat_distance(lon1, lat1, lon2, lat2);
    }
  }
  return mat;
}

} // namespace EpiGraph
#endif // EPIGRAPH_SPATIAL_UTIL_H
