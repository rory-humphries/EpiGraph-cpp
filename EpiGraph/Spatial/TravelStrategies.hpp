//
// Created by roryh on 14/10/2020.
//

#ifndef EPIGRAPH_TRAVEL_STRATEGIES_H
#define EPIGRAPH_TRAVEL_STRATEGIES_H

#include <Eigen/Core>
namespace EpiGraph {
template <typename Mat1, typename Mat2, typename Mat3>
auto update_out_travel_weights_with_compliance(
    Eigen::MatrixBase<Mat1> &travel_weights, Eigen::ArrayBase<Mat2> &bool_mat,
    Eigen::MatrixBase<Mat3> &compliance_vec) -> Mat1 {

  Mat1 new_travel_weights = (bool_mat).select(
      travel_weights.array(),
      travel_weights.array() *
          (1 -
           compliance_vec.rowwise().replicate(travel_weights.cols()).array()));
  return new_travel_weights;
}
} // namespace EpiGraph
#endif // EPIGRAPH_TRAVEL_STRATEGIES_H
