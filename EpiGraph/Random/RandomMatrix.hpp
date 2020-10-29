//
// Created by roryh on 24/06/2020.
//

#ifndef EPIGRAPH_RANDOM_MATRIX_H
#define EPIGRAPH_RANDOM_MATRIX_H

#include <vector>
#include <random>
#include <memory>
#include <algorithm>
#include <Eigen/Sparse>
#include <Eigen/Dense>

#include <EpiGraph/Random/Distributions.hpp>
#include <EpiGraph/Core/IO.hpp>

namespace EpiGraph {

    class RandomMatrixGenerator {
    private:
        size_t dim;
        std::vector<std::discrete_distribution<int>> m_row_distributions;

    public:

        explicit RandomMatrixGenerator(size_t dim) : dim(dim), m_row_distributions(dim) {}


        template<typename TIter>
        auto set_row_distribution(size_t i, TIter begin, TIter end) -> void {
            // stores the probability distributions of a traveller going to a destination vertex given the traveller is
            // leaving from vertex v_src. The index of the vector correspond to the vertices of the same index. i.e. vec[v_src]

            if (std::distance(begin, end) != dim)
                throw std::invalid_argument("Iterator length does not match dimension of model");
            if (i > dim)
                throw std::invalid_argument("Index greater than dimension of model");

            m_row_distributions[i] = std::discrete_distribution<int>(begin, end);
        }

        template<typename Derived>
        auto set_row_distributions(const Eigen::DenseBase<Derived> &weights) -> void {
            if (weights.cols() != dim)
                throw std::invalid_argument("Matrix cols do not match model dimension");

#pragma omp parallel for
            for (int i = 0; i < weights.rows(); i++) {
                Eigen::RowVectorXd row = weights.row(i);
                set_row_distribution(i, row.data(), row.data() + weights.cols());
            }
        }

        template<typename DerivedB>
        auto
        distribute_vec_over_matrix_rows(Eigen::MatrixBase<DerivedB> &vals) -> Eigen::SparseMatrix<double> {
            /*
             * Add random edges to the network with a random number of travellers along each edge. The total number of
             * travelers out of v_src is decided by the commuter_dist distribution which gives the proportion of the
             * total population in v_src that will travel.
             *
             * travel distribution is a random distribution that supports the operator() which take a random generator and
             * outputs a destination vertex.
             */

            Eigen::SparseMatrix<double> adj(vals.size(), vals.size());

            std::vector<Eigen::Triplet<double>> tripletList;
            tripletList.reserve(vals.sum());
            std::random_device rd;
            std::vector<std::mt19937> gen_vec;
            for (int i = 0; i < omp_get_max_threads(); i++) gen_vec.emplace_back(rd());

#pragma omp parallel for
            for (int i = 0; i < dim; i++) {
                std::discrete_distribution<int> &travel_distribution = m_row_distributions[i];
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
}

#endif //EPIGRAPH_RANDOM_MATRIX_H
