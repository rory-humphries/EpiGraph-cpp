
#ifndef EPIGRAPH_BFS_HPP
#define EPIGRAPH_BFS_HPP

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <EpiGraph/EigenUtil/StaticAsserts.hpp>
#include <iostream>
#include <unordered_set>

namespace EpiGraph {

template <IsSparseOrDenseMatrix Adj>
auto BFS(const Adj &adj, Eigen::Index v)
    -> int /*Eigen::SparseMatrix<typename Adj::Scalar>*/ {
  Eigen::SparseVector<typename Adj::Scalar> x(adj.cols());
  Eigen::SparseVector<typename Adj::Scalar> xold(adj.cols());

  std::unordered_set<int> found;

  Adj adj_cp = adj;

  x.coeffRef(v) = 1;
  int nnz = 1;
  while (x.nonZeros() != xold.nonZeros()) {

    nnz = x.nonZeros();
    // std::cout << "x: " << x.nonZeros() << "xold: " << xold.nonZeros()
    //          << std::endl;
    xold = x;
    x = x + (adj * x);
    // std::cout << "x: " << x.nonZeros() << "xold: " << xold.nonZeros()
    //          << std::endl;

    /*
        for (typename Eigen::SparseVector<typename Adj::Scalar>::InnerIterator
       it( xold); it; ++it) {

          std::cout << it.index() << std::endl;
          adj_cp.col(it.index()) *= 0;
          x.coeffRef(it.index()) = 0;
          found.insert(it.index());
        }
        */
  }
  return x.nonZeros();
}

template <IsSparseMatrix Adj>
auto BFS2(const Adj &adj, Eigen::Index v)
    -> int /*Eigen::SparseMatrix<typename Adj::Scalar>*/ {
  Eigen::VectorXi T = Eigen::VectorXi::Zero(adj.cols());
  Eigen::VectorXi L = Eigen::VectorXi::Zero(adj.cols());

  int start = 0;
  int end = 1;
  int z = 1;

  L[start] = v;

  while (start != end) {
    for (int j = start; j < end; j++) {
      int node = L[j];
      T[node] = 1;
      for (typename Adj::InnerIterator it(adj, node); it; ++it) {
        int i = it.row();

        if (T[i] == 0) {
          T[i] = 1;
          L[z] = i;
          z++;
        }
      }
    }
    start = end;
    end = z;
  }
  return T.sum();
}

template <IsSparseMatrix Adj>
auto BFS_parallel(const Adj &adj, Eigen::Index v)
    -> int /*Eigen::SparseMatrix<typename Adj::Scalar>*/ {
  Eigen::VectorXi T = Eigen::VectorXi::Zero(adj.cols());

  Eigen::VectorXi x = Eigen::VectorXi::Zero(adj.cols());

  std::vector<int> L;
  L.reserve(adj.cols());

  int start = 0;
  int end = 1;
  int z = 1;

  L[start] = v;
  // x[v] = 1;
  // T[v] = 1;

  std::vector<std::vector<int>> pvec(16);

  while (start != end) {
#pragma omp parallel for
    for (int j = start; j < end; j++) {

      int node = L[j];
      T[node] = 1;

      int tid = omp_get_thread_num();

      for (typename Adj::InnerIterator it(adj, node); it; ++it) {
        int i = it.row();

        // x[ind] = 1;
        if (T[i] == 0) {
          pvec[tid].push_back(i);
        }
      }
    }

    for (auto &v : pvec) {
      for (auto &i : v) {
        if (T[i] == 0) {
          T[i] = 1;
          L[z] = i;
          z++;
        }
        v.clear();
      }
    }
    start = end;
    end = z;
  }
  return T.sum();
}

} // namespace EpiGraph

#endif
