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

class RandomTravelModel {
public:
    double m_compliance;
    double m_max_dist;

    std::vector<std::vector<double>> rand_mat;
    ProbDist commuter_dist;

    RandomTravelModel(std::string rand_mat_path, std::string commuter_dist_path) :
            rand_mat(matrix_from_csv(rand_mat_path)) {
        commuter_dist = ProbDist_from_csv(commuter_dist_path);
    }


    auto add_random_edges(SpatialMetaPopNetwork &g) -> void {

        std::uniform_real_distribution<> uni_dist(0, 1);

        int travs = 0;
        double u;
        for (int v_src = 0; v_src < g.graph.num_vertices(); v_src++) {

            // calculate number of commuters
            double travel_prop = commuter_dist(global_engine());
            int N = g.populations[v_src] * travel_prop;

            std::map<Vertex, double> edges_to_add;
            for (int n = 0; n < N; n++) {
                // find the destination vertex from the random matrix
                u = uni_dist(global_engine());
                auto it = std::lower_bound(rand_mat[v_src].begin(), rand_mat[v_src].end(), u);
                Vertex v_dst = std::distance(rand_mat[v_src].begin(), it);

                // if probabilities don't exactly add to one, v_dst can be one larger than the total number of vertices
                v_dst = std::min((int) v_dst, (int) g.get_graph().num_vertices() - 1);

                if (v_dst == v_src)
                    continue;

                // if distance is greater than max distance only non compliant will travel
                if (long_lat_distance(g.longlat[v_src], g.longlat[v_dst]) >= m_max_dist) {
                    u = uni_dist(global_engine());
                    if (u >= m_compliance)
                        continue;
                }
                edges_to_add[v_dst] += 1;
            }
            for (auto &k: edges_to_add) {
                g.add_coupling(v_src, k.first, k.second);
            }
        }
    }
};

#endif //EPIGRAPH_TRAVEL_MODEL_H
