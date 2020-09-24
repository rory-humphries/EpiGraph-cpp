//
// Created by roryh on 24/06/2020.
//

#ifndef EPIGRAPH_TRAVEL_MODEL_H
#define EPIGRAPH_TRAVEL_MODEL_H

#include <vector>
#include <random>
#include <memory>
#include <algorithm>
#include <Eigen/Sparse>
#include <Eigen/Dense>

#include "distributions.h"
#include "io.h"

class RandomMatrix {
private:
    size_t dim;
    std::vector<std::discrete_distribution<int>> m_entry_probs;

public:

    RandomMatrix(size_t dim) : dim(dim), m_entry_probs(dim) {}


    template<typename TIter>
    auto set_entry_probs(size_t i, TIter begin, TIter end) -> void {
        // stores the probability distributions of a traveller going to a destination vertex given the traveller is
        // leaving from vertex v_src. The index of the vector correspond to the vertices of the same index. i.e. vec[v_src]

        if (std::distance(begin, end) != dim)
            throw std::invalid_argument("Iterator length does not match dimension of model");
        if (i > dim)
            throw std::invalid_argument("Index greater than dimension of model");

        m_entry_probs[i] = std::discrete_distribution<int>(begin, end);
    }

    template<typename Derived>
    auto set_entry_probs(Eigen::DenseBase<Derived> &weights) -> void {
        if (weights.cols() != dim)
            throw std::invalid_argument("Matrix cols do not match model dimension");

        for (int i = 0; i < weights.rows(); i++) {
            Eigen::RowVectorXd row = weights.row(i);
            set_entry_probs(i, row.data(), row.data() + weights.cols());
        }
    }

    template<typename DerivedA, typename DerivedB>
    auto generate(Eigen::SparseMatrix<DerivedA> &adj, Eigen::MatrixBase<DerivedB> &vals) -> void {
        /*
         * Add random edges to the network with a random number of travellers along each edge. The total number of
         * travelers out of v_src is decided by the commuter_dist distribution which gives the proportion of the
         * total population in v_src that will travel.
         *
         * travel distribution is a random distribution that supports the operator() which take a random generator and
         * outputs a destination vertex.
         */

        //Eigen::SparseMatrix<double> adj_tmp(pop_map.size(), pop_map.size());
        std::uniform_real_distribution<> uni_dist(0, 1);

        std::vector<Eigen::Triplet<double>> tripletList;
        tripletList.reserve(vals.sum());
        std::map<std::pair<int, int>, double> edges_to_add;

        for (int i = 0; i < dim; i++) {
            std::discrete_distribution<int> &travel_distribution = m_entry_probs[i];

            for (int n = 0; n < vals[i]; n++) {
                // find the destination vertex from the travel distribution
                int j = travel_distribution(global_engine());

                tripletList.emplace_back(i, j, 1);
            }
        }

        adj.resize(dim, dim);
        adj.setFromTriplets(tripletList.begin(), tripletList.end());
        adj.makeCompressed();
    }

};


#endif //EPIGRAPH_TRAVEL_MODEL_H
