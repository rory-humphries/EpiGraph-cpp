//
// Created by roryh on 24/06/2020.
//

#ifndef EPIGRAPH_RANDOM_MATRIX_H
#define EPIGRAPH_RANDOM_MATRIX_H

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <algorithm>
#include <memory>
#include <random>
#include <vector>

namespace EpiGraph {

/**
 * @brief A container for managing a list of std::discrete_distributions<int>
 * which all produce random integers in the same range with different
 * probabilities.
 *
 */
class RandomAdjMat {
private:
  size_t dim;
  std::vector<std::discrete_distribution<int>> m_distributions_vector;

public:
  RandomAdjMat() : dim(0), m_distributions_vector(0) {}

  explicit RandomAdjMat(size_t dim) : dim(dim), m_distributions_vector(dim) {}

  template <typename TIter>
  auto set_distribution(size_t i, TIter begin, TIter end) -> void {
    if (std::distance(begin, end) != dim)
      throw std::invalid_argument(
          "Iterator length does not match dimension of model");
    if (i > dim)
      throw std::invalid_argument("Index greater than dimension of model");

    m_distributions_vector[i] = std::discrete_distribution<int>(begin, end);
  }

  template <typename Derived>
  auto set_distributions(const Eigen::DenseBase<Derived> &weights) -> void {
    if (weights.cols() != dim)
      throw std::invalid_argument("Matrix cols do not match model dimension");

#pragma omp parallel for
    for (int i = 0; i < weights.rows(); i++) {
      Eigen::RowVectorXd row = weights.row(i);
      set_distribution(i, row.data(), row.data() + weights.cols());
    }
  }

  auto get_vector() const
      -> const std::vector<std::discrete_distribution<int>> & {
    return m_distributions_vector;
  }

  template <typename Derived>
  auto gen_sparse_mat(Eigen::MatrixBase<Derived> &vals)
      -> Eigen::SparseMatrix<double> {
    /*
     * Add random edges to the network with a random number of travellers along
     * each edge. The total number of travelers out of v_src is decided by the
     * commuter_dist distribution which gives the proportion of the total
     * population in v_src that will travel.
     *
     * travel distribution is a random distribution that supports the operator()
     * which take a random generator and outputs a destination vertex.
     */

    Eigen::SparseMatrix<double> adj(vals.size(), vals.size());

    std::vector<Eigen::Triplet<double>> tripletList;
    tripletList.reserve(vals.sum());
    std::random_device rd;
    std::vector<std::mt19937> gen_vec;
    for (int i = 0; i < omp_get_max_threads(); i++)
      gen_vec.emplace_back(rd());

#pragma omp parallel for
    for (int i = 0; i < dim; i++) {
      std::discrete_distribution<int> &travel_distribution =
          m_distributions_vector[i];
      std::mt19937 &gen = gen_vec[omp_get_thread_num()];

      std::vector<Eigen::Triplet<double>> p_vec;
      p_vec.reserve(vals[i]);
      for (int n = 0; n < vals[i]; n++) {

        // find the destination vertex from the travel distribution
        int j = travel_distribution(gen);
        p_vec.emplace_back(i, j, 1);
      }
#pragma omp critical
      tripletList.insert(tripletList.end(), p_vec.begin(), p_vec.end());
    }

    adj.setFromTriplets(tripletList.begin(), tripletList.end());
    adj.makeCompressed();
    return adj;
  }
};

} // namespace EpiGraph

#endif // EPIGRAPH_RANDOM_MATRIX_H
