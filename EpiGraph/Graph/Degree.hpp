
#ifndef EPIGRAPH_DEGREE_HPP
#define EPIGRAPH_DEGREE_HPP

#include <Eigen/Core>
#include <EpiGraph/EigenUtil/StaticAsserts.hpp>

namespace EpiGraph {

/*
 * Undirected degree
 */
template <IsSparseOrDenseMatrix Adj>
auto deg(const Adj &adj, Eigen::Index node) -> typename Adj::Scalar {

  return adj.block(node, node, 1, adj.cols() - node).sum() +
         adj.block(0, node, node, 1).sum();
  //-------------------------^ node +1 will double count self loops
}

template <IsSparseOrDenseMatrix Adj>
auto deg(const Adj &adj)
    -> Eigen::Matrix<typename Adj::Scalar, Adj::RowsAtCompileTime, 1> {
  using Vec =
      typename Eigen::Matrix<typename Adj::Scalar, Adj::RowsAtCompileTime, 1>;

  Eigen::TriangularView viewU = adj.template triangularView<Eigen::Upper>();
  Eigen::TriangularView viewSU =
      adj.template triangularView<Eigen::StrictlyUpper>();

  return (Vec::Ones(adj.rows(), 1).transpose() * viewU) +
         (viewSU * Vec::Ones(adj.rows(), 1)).transpose();
}

/*
 * Out degree
 */

template <IsMatrix Adj>
auto out_deg(const Adj &adj)
    -> Eigen::Matrix<typename Adj::Scalar, Adj::RowsAtCompileTime, 1> {
  return adj.colwise().sum().transpose();
}

template <IsSparseMatrix Adj>
auto out_deg(const Adj &adj)
    -> Eigen::Matrix<typename Adj::Scalar, Adj::RowsAtCompileTime, 1> {
  return (Eigen::Matrix<typename Adj::Scalar, 1, Adj::ColsAtCompileTime>::Ones(
              1, adj.cols()) *
          adj.derived())
      .transpose();
}

/*
 * In degree
 */

template <IsMatrix Adj>
auto in_deg(const Adj &adj)
    -> Eigen::Matrix<typename Adj::Scalar, Adj::RowsAtCompileTime, 1> {
  return adj.rowwise().sum();
}

template <IsSparseMatrix Adj>
auto in_deg(const Adj &adj)
    -> Eigen::Matrix<typename Adj::Scalar, Adj::RowsAtCompileTime, 1> {
  using Vec =
      typename Eigen::Matrix<typename Adj::Scalar, Adj::RowsAtCompileTime, 1>;

  return adj.derived() * Vec::Ones(adj.rows(), 1);
}

} // namespace EpiGraph

#endif
