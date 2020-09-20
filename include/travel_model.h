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
#include "spatial_util.h"

class RandomTravelModel {
private:
    size_t dim;

    double m_compliance;
    double m_max_dist;

    std::vector<std::discrete_distribution<int>> travel_probs;

public:
    ProbDist commuter_dist;
    ProbDist travel_dist;

public:

    RandomTravelModel(size_t dim) : dim(dim), travel_probs(dim), commuter_dist(), travel_dist(), m_compliance(0),
                                    m_max_dist(0) {}

    // getters
    auto get_compliance(double compliance) -> double {
        return m_compliance;
    }

    auto get_max_dist(double max_dist) -> double {
        return m_max_dist;
    }

    // setters
    auto set_compliance(double compliance) -> void {
        m_compliance = compliance;
    }

    auto set_max_distance(double max_dist) -> void {
        m_max_dist = max_dist;
    }

    auto set_commuter_dist(ProbDist commuter_dist) -> void {
        this->commuter_dist = commuter_dist;
    }

    template<typename TIter>
    auto set_travel_probs(size_t i, TIter begin, TIter end) -> void {
        // stores the probability distributions of a traveller going to a destination vertex given the traveller is
        // leaving from vertex v_src. The index of the vector correspond to the vertices of the same index. i.e. vec[v_src]

        if (std::distance(begin, end) != dim)
            throw std::invalid_argument("Iterator length does not match dimension of model");
        if (i > dim)
            throw std::invalid_argument("Index greater than dimension of model");

        travel_probs[i] = std::discrete_distribution<int>(begin, end);
    }
    auto set_travel_probs(std::string path_to_csv) {
        size_t i = 0;
        std::ifstream infile(path_to_csv);
        std::string line;
        while (getline(infile, line, '\n')) {
            std::vector<double> tmp;
            std::istringstream ss(line);
            std::string token;

            while (std::getline(ss, token, ',')) {
                tmp.push_back(stof(token));
            }
            set_travel_probs(i, tmp.begin(), tmp.end());
            i++;
        }
    }


    template<typename Derived, typename DerivedB, typename DerivedC>
    auto add_out_travels(Eigen::SparseMatrix<Derived> &adj, Eigen::MatrixBase<DerivedB> &pos_map,
                         Eigen::MatrixBase<DerivedC> &pop_map) -> void {
        /*
         * Add random edges to the network with a random number of travellers along each edge. The total number of
         * travelers out of v_src is decided by the commuter_dist distribution which gives the proportion of the
         * total population in v_src that will travel.
         *
         * travel distribution is a random distribution that supports the operator() which take a random generator and
         * outputs a destination vertex.
         */

        int tot_pop = pop_map.sum();

        //Eigen::SparseMatrix<double> adj_tmp(pop_map.size(), pop_map.size());
        std::uniform_real_distribution<> uni_dist(0, 1);

        typedef Eigen::Triplet<double> T;
        std::vector<T> tripletList;
        tripletList.reserve(tot_pop);
        std::map<std::pair<int, int>, double> edges_to_add;

        for (int i = 0; i < dim; i++) {
            std::discrete_distribution<int> &travel_distribution = travel_probs[i];

            // get the number of people travelling
            double travel_prop = commuter_dist(global_engine());
            int N = pop_map[i] * travel_prop;

            for (int n = 0; n < N; n++) {

                // find the destination vertex from the travel distribution
                int j = travel_distribution(global_engine());

                // if distance is greater than max distance only non compliant will travel
                double lon1 = pos_map(i, 0), lon2 = pos_map(j, 0), lat1 = pos_map(i, 1), lat2 = pos_map(j, 1);

                if (long_lat_distance(lon1, lat1, lon2, lat2) / 1000.0 >= m_max_dist) {
                    double u = uni_dist(global_engine());
                    if (u >= m_compliance) {
                        tripletList.emplace_back(i, j, 1);
                    }
                } else {
                    tripletList.emplace_back(i, j, 1);
                }
            }
        }

        adj.setFromTriplets(tripletList.begin(), tripletList.end());
        adj.makeCompressed();
    }

};


#endif //EPIGRAPH_TRAVEL_MODEL_H
