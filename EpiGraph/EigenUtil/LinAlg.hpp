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
auto general_matrix_prod(EigenUtil::MatrixBase<DerivedA> &A, EigenUtil::MatrixBase<DerivedB> &B,
                         std::function<typename DerivedA::Scalar (EigenUtil::MatrixBase<DerivedA> &A_row,
                                                                 EigenUtil::MatrixBase<DerivedB> &B_col)> prod) -> void {
    col_vector_assert(vec);
    col_vector_assert(groups);

    op_vec = EigenUtil::Matrix<typename DerivedA::Scalar, EigenUtil::Dynamic, 1>::Zero(groups.maxCoeff() + 1);

    for (int i = 0; i < vec.rows(); i++) {
        op_vec[groups[i]] += vec[i];
    }
}
*/
}
#endif //EPIGRAPH_EIGEN_UTIL_H
