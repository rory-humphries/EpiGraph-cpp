//
// Created by roryh on 19/09/2020.
//

#ifndef EPIGRAPH_EIGEN_UTIL_H
#define EPIGRAPH_EIGEN_UTIL_H

#include <Eigen/Core>
#include <fstream>
#include "spatial_util.h"

template<typename Derived>
auto col_vector_assert(const Eigen::MatrixBase<Derived> &vec) -> void {
    static_assert(Derived::ColsAtCompileTime == 1, "Expected a column vector");
}

template<typename Derived>
auto row_vector_assert(const Eigen::MatrixBase<Derived> &vec) -> void {
    static_assert(Derived::RowsAtCompileTime == 1, "Expected a row vector");
}

template<typename DerivedA, typename DerivedB>
auto distance_matrix(Eigen::MatrixBase<DerivedB> &pos_mat) -> DerivedA {
    static_assert(DerivedB::ColsAtCompileTime == 2, "Expected a N x 2 matrix");

    DerivedA mat(pos_mat.rows(), pos_mat.rows());
    mat.setZero();

    for (int i = 0; i < pos_mat.rows(); i++) {
        for (int j = 0; j < pos_mat.rows(); j++) {

            double lon1 = pos_mat(i, 0), lon2 = pos_mat(j, 0), lat1 = pos_mat(i, 1), lat2 = pos_mat(j, 1);
            mat(i, j) = long_lat_distance(lon1, lat1, lon2, lat2);
        }
    }
    return mat;
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

template<typename Derived, typename T>
auto accumulate_groups(const Eigen::MatrixBase<Derived> &vec1,
                       const std::vector<T> &vec2) -> std::map<T, typename Derived::Scalar> {
    /*
     * For each unique value in vec2, accumulate the elements in vec1 which match indices with all the occurrences of a
     * given element in vec2.
     * TODO: Add error checks
     */
    std::map<T, typename Derived::Scalar> group_sums;
    for (int i = 0; i < vec1.rows(); i++) group_sums[vec2[i]] += vec1[i];

    return group_sums;
}

template<typename DerivedA, typename DerivedB>
auto mat_rowwise_less_than_vec(const Eigen::ArrayBase<DerivedA> &mat,
                               const Eigen::ArrayBase<DerivedB> &vec) -> Eigen::Array<
        bool,
        DerivedA::RowsAtCompileTime,
        DerivedA::ColsAtCompileTime> {
    return mat < vec.rowwise().replicate(mat.cols());
}

template<typename DerivedA, typename DerivedB>
auto mat_colwise_less_than_vec(const Eigen::ArrayBase<DerivedA> &mat,
                               const Eigen::ArrayBase<DerivedB> &vec) -> Eigen::Array<
        bool,
        DerivedA::RowsAtCompileTime,
        DerivedA::ColsAtCompileTime> {
    return mat < vec.colwise().replicate(mat.rows());
}

template<typename Derived, typename F>
auto general_outer_product(
        const Eigen::MatrixBase<Derived> &v,
        const Eigen::MatrixBase<Derived> &u,
        F f) -> Eigen::Matrix<typename Derived::Scalar, Derived::RowsAtCompileTime,
        Derived::RowsAtCompileTime> {
    col_vector_assert(v);
    col_vector_assert(u);
    Eigen::Matrix<typename Derived::Scalar, Derived::RowsAtCompileTime, Derived::RowsAtCompileTime> mat(v.rows(),
                                                                                                        u.rows());
    for (int j = 0; j < v.rows(); j++) {
        for (int i = 0; i < u.rows(); i++) {
            mat(i, j) = f(v[i], u[j]);
        }
    }
    return mat;
}

/*
template<typename DerivedA, typename DerivedB>
auto general_matrix_prod(Eigen::MatrixBase<DerivedA> &A, Eigen::MatrixBase<DerivedB> &B,
                         std::function<typename DerivedA::Scalar (Eigen::MatrixBase<DerivedA> &A_row,
                                                                 Eigen::MatrixBase<DerivedB> &B_col)> prod) -> void {
    col_vector_assert(vec);
    col_vector_assert(groups);

    op_vec = Eigen::Matrix<typename DerivedA::Scalar, Eigen::Dynamic, 1>::Zero(groups.maxCoeff() + 1);

    for (int i = 0; i < vec.rows(); i++) {
        op_vec[groups[i]] += vec[i];
    }
}
*/
#endif //EPIGRAPH_EIGEN_UTIL_H
