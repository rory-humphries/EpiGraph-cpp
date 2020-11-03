//
// Created by roryh on 02/11/2020.
//

#ifndef EPIGRAPH_CPP_REDUCTIONS_HPP
#define EPIGRAPH_CPP_REDUCTIONS_HPP

#include <EpiGraph/EigenUtil/StaticAsserts.hpp>

#include <Eigen/Dense>

#include <map>
#include <vector>


namespace EpiGraph {
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
}
#endif //EPIGRAPH_CPP_REDUCTIONS_HPP
