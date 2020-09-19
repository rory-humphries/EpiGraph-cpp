//
// Created by roryh on 19/09/2020.
//

#ifndef EPIGRAPH_EIGEN_UTIL_H
#define EPIGRAPH_EIGEN_UTIL_H

#include <Eigen/Core>
#include <fstream>

template<typename Derived>
auto col_vector_assert(Eigen::MatrixBase<Derived>& vec) -> void {
    static_assert(Derived::ColsAtCompileTime == 1, "Expected a column vector");
}

template<typename Derived>
auto row_vector_assert(Eigen::MatrixBase<Derived>& vec) -> void {
    static_assert(Derived::RowsAtCompileTime == 1, "Expected a row vector");
}


#endif //EPIGRAPH_EIGEN_UTIL_H
