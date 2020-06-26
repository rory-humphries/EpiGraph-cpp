//
// Created by roryh on 24/06/2020.
//

#ifndef EPIGRAPH_TRAVEL_MODEL_H
#define EPIGRAPH_TRAVEL_MODEL_H

#include <vector>
#include <random>
#include "distributions.h"
#include "meta_pop_network.h"
#include "io.h"
#include "spatial_utils.h"

template<typename TNetwork>
class RandomTravelModel {
public:
    double m_compliance;
    double m_max_dist;

    TNetwork* network;

    ProbDist commuter_dist;
    std::vector<std::discrete_distribution<int>> travel_dist;
    RandomTravelModel(std::string rand_mat_path, std::string commuter_dist_path) :
            commuter_dist(ProbDist_from_csv(commuter_dist_path)) {

        auto vec = matrix_from_csv(rand_mat_path);
        for (auto &k: vec) {
            travel_dist.emplace_back(k.begin(), k.end());
        }
    }

    auto set_network(TNetwork& g) -> void{
        network = &g;
    }

    auto add_random_edges() -> void {

        std::uniform_real_distribution<> uni_dist(0, 1);

        double u;
        for (int v_src = 0; v_src < (*network).num_vertices(); v_src++) {

            // calculate number of commuters
            double travel_prop = commuter_dist(global_engine());
            int N = (*network).vprop[v_src].population * travel_prop;

            std::map<Vertex, double> edges_to_add;
            for (int n = 0; n < N; n++) {

                // find the destination vertex from the random matrix
                u = uni_dist(global_engine());
                Vertex v_dst = travel_dist[v_src](global_engine());
                if (v_dst == v_src)
                    continue;

                // if distance is greater than max distance only non compliant will travel
                if (long_lat_distance((*network).vprop[v_src].position,
                        (*network).vprop[v_dst].position) >= m_max_dist) {
                    u = uni_dist(global_engine());
                    if (u >= m_compliance)
                        continue;
                }
                edges_to_add[v_dst] += 1;
            }
            for (auto &k: edges_to_add) {
                // maybe problems if network is ever undirected
                auto e_pair = (*network).add_edge(v_src, k.first);
                Edge e = e_pair.first;

                (*network).eprop[e].population = k.second;
            }
        }
    }
};

#endif //EPIGRAPH_TRAVEL_MODEL_H
