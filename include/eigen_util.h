//
// Created by roryh on 19/09/2020.
//

#ifndef EPIGRAPH_EIGEN_UTIL_H
#define EPIGRAPH_EIGEN_UTIL_H

#include <Eigen/Dense>
#include <fstream>
#include "spatial_util.h"

template<typename Derived>
auto col_vector_assert(Eigen::MatrixBase<Derived> &vec) -> void {
    static_assert(Derived::ColsAtCompileTime == 1, "Expected a column vector");
}

template<typename Derived>
auto row_vector_assert(Eigen::MatrixBase<Derived> &vec) -> void {
    static_assert(Derived::RowsAtCompileTime == 1, "Expected a row vector");
}

template<typename DerivedA, typename DerivedB>
auto distance_matrix(Eigen::MatrixBase<DerivedA> &mat, Eigen::MatrixBase<DerivedB> &pos_mat) -> void {
    static_assert(DerivedB::ColsAtCompileTime == 2, "Expected a N x 2 matrix");

    mat = DerivedA(pos_mat.rows(), pos_mat.rows());
    mat.setZero();

    for (int i = 0; i < pos_mat.rows(); i++) {
        for (int j = 0; j < pos_mat.rows(); j++) {

            double lon1 = pos_mat(i, 0), lon2 = pos_mat(j, 0), lat1 = pos_mat(i, 1), lat2 = pos_mat(j, 1);
            mat(i, j) = long_lat_distance(lon1, lat1, lon2, lat2) / 1000.0;
        }
    }
}

template<typename DerivedA, typename DerivedB, typename DerivedC>
auto sum_groups(Eigen::MatrixBase<DerivedA> &op_vec, Eigen::MatrixBase<DerivedB> &vec,
                Eigen::MatrixBase<DerivedC> &groups) -> void {
    col_vector_assert(vec);
    col_vector_assert(groups);

    op_vec = Eigen::Matrix<typename DerivedA::Scalar, Eigen::Dynamic, 1>::Zero(groups.maxCoeff() + 1);

    for (int i = 0; i < vec.rows(); i++) {
        op_vec[groups[i]] += vec[i];
    }
}

#endif //EPIGRAPH_EIGEN_UTIL_H
