//
// Created by roryh on 14/10/2020.
//

#ifndef EPIGRAPH_TRAVEL_STRATEGIES_H
#define EPIGRAPH_TRAVEL_STRATEGIES_H

#include <Eigen/Core>

template<typename DerivedA, typename DerivedB, typename DerivedC>
auto update_out_travel_weights_with_compliance(Eigen::MatrixBase<DerivedA> &travel_weights,
                                               Eigen::ArrayBase<DerivedB> &bool_mat,
                                               Eigen::MatrixBase<DerivedC> &compliance_vec) -> DerivedA {

    DerivedA new_travel_weights = (bool_mat).select(travel_weights.array(),
                                                    travel_weights.array() *
                                                    (1 - compliance_vec.rowwise().replicate(
                                                            travel_weights.cols()).array()));
    return new_travel_weights;
}

#endif //EPIGRAPH_TRAVEL_STRATEGIES_H
