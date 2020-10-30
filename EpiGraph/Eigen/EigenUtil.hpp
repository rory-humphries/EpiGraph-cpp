//
// Created by roryh on 19/09/2020.
//

#ifndef EPIGRAPH_EIGEN_UTIL_H
#define EPIGRAPH_EIGEN_UTIL_H

#include <Eigen/Core>
#include <fstream>
#include <map>
#include <EpiGraph/Spatial/SpatialUtil.hpp>

namespace EpiGraph {
    template<typename Derived>
    auto vector_assert(const Eigen::MatrixBase<Derived> &vec) -> void {
        static_assert(Derived::IsVectorAtCompileTime, "Expected a vector");
    }

    template<typename Derived>
    auto col_vector_assert(const Eigen::MatrixBase<Derived> &vec) -> void {
        static_assert(Derived::ColsAtCompileTime == 1, "Expected a column vector");
    }

    template<typename Derived>
    auto row_vector_assert(const Eigen::MatrixBase<Derived> &vec) -> void {
        static_assert(Derived::RowsAtCompileTime == 1, "Expected a row vector");
    }

    template<typename Derived, typename T>
    auto accumulate_groups(const Eigen::MatrixBase<Derived> &vec1,
                           const std::vector<T> &vec2) -> std::map<T, typename Derived::Scalar> {
        /*
         * For each unique value in vec2, accumulate the elements in vec1 which match indices with all the occurrences of a
         * given element in vec2.
         */
        col_vector_assert(vec1);

        std::map<T, typename Derived::Scalar> group_sums;
        for (int i = 0; i < vec1.rows(); i++) group_sums[vec2[i]] += vec1[i];

        return group_sums;
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
}
#endif //EPIGRAPH_EIGEN_UTIL_H
