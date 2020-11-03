//
// Created by roryh on 03/11/2020.
//

#ifndef EPIGRAPH_CPP_STATICASSERTS_HPP
#define EPIGRAPH_CPP_STATICASSERTS_HPP

#include <Eigen/Dense>

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
}
#endif //EPIGRAPH_CPP_STATICASSERTS_HPP
