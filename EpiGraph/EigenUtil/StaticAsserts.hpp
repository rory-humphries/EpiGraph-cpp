//
// Created by roryh on 03/11/2020.
//

#ifndef EPIGRAPH_CPP_STATICASSERTS_HPP
#define EPIGRAPH_CPP_STATICASSERTS_HPP

#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace EpiGraph {

template <typename Derived>
auto vector_assert(const Eigen::MatrixBase<Derived> &vec) -> void {
  static_assert(Derived::IsVectorAtCompileTime, "Expected a vector");
}

template <typename Derived>
auto col_vector_assert(const Eigen::MatrixBase<Derived> &vec) -> void {
  static_assert(Derived::ColsAtCompileTime == 1, "Expected a column vector");
}

template <typename Derived>
auto row_vector_assert(const Eigen::MatrixBase<Derived> &vec) -> void {
  static_assert(Derived::RowsAtCompileTime == 1, "Expected a row vector");
}

template <typename T>
concept IsMatrix = std::is_base_of<Eigen::MatrixBase<T>, T>::value;

template <typename T>
concept IsSparseMatrix = std::is_base_of<Eigen::SparseMatrixBase<T>, T>::value;

template <typename T>
concept IsSparseOrDenseMatrix = IsMatrix<T> || IsSparseMatrix<T>;

template <typename T>
concept IsColVector = IsMatrix<T> &&T::ColsAtCompileTime == 1;

template <typename T>
concept IsRowVector = IsMatrix<T> &&T::RowsAtCompileTime == 1;

template <typename T> concept IsVector = IsColVector<T> || IsRowVector<T>;

} // namespace EpiGraph
#endif // EPIGRAPH_CPP_STATICASSERTS_HPP
