//
// Created by roryh on 23/06/2020.
//

#ifndef EPIGRAPH_META_POP_NETWORK_H
#define EPIGRAPH_META_POP_NETWORK_H

#include <vector>
#include <map>
#include "edge.h"
#include "graph.h"

struct VMetaPop {
    /*
     * A vertex property that holds a population
     */
    using population_type = double;

    population_type population;
};

struct EMetaPop {
    /*
     * An edge property that holds a population
     */
    using population_type = double;

    population_type population;
};

struct VSpatialMetaPop {
    /*
     * A vertex property that holds a population and location
     */
    using population_type = double;
    using position_type = std::pair<double, double>;

    population_type population;
    position_type position;
};

class SpatialMetaPopNetwork {

public:
    using coupling_type = std::map<Edge, double>;
    using OutEdgeIter = Graph::OutEdgeIter;
    using InEdgeIter = Graph::OutEdgeIter;
    using EdgeIter = Graph::EdgeIter;

    // the underlying graph structure
    Graph graph;

    // holds the population of each meta population
    std::vector<double> populations;

    // holds the longlat position of each meta population
    std::vector<std::pair<double, double>> longlat;

    // holds the number of travellers between each metapopulation
    coupling_type weights;

    // constructors
    SpatialMetaPopNetwork() : graph() {}


    auto get_graph() const -> const Graph & {
        /*
         * Return a const reference to the underlying graph structure.
         */
        return graph;
    }

    auto add_metapop() -> Vertex {
        /*
         * Add a metapopulation to the network with population 0 located at 0,0.
         */
        auto v = graph.add_vertex();
        populations.emplace_back();
        longlat.emplace_back();

        return v;
    }

    auto add_metapop(double N, double longatude, double latatude) -> void {
        /*
         * Add a metapopulation to the network with population N located at longatude, latatude.
         */
        graph.add_vertex();
        populations.push_back(N);
        longlat.emplace_back(longatude, latatude);
    }

    auto out_edges(Vertex v) -> std::pair<OutEdgeIter, OutEdgeIter> {
        return graph.out_edges(v);
    }

    auto set_population(Vertex v, double N) -> void {
        populations.at(v) = N;
    }
    auto get_population(Vertex v) -> double {
        return populations.at(v);
    }
    auto set_longlat(Vertex v, double lo, double la) {
        longlat.at(v) = {lo, la};
    }
    auto get_longlat(Vertex v) -> std::pair<double, double> {
        return longlat.at(v);
    }

    auto add_coupling(int vsrc, int vdst, double weight, bool sum_weights = false) -> void {
        /*
         * If the edge between vsrc and vdst does not exist, create it and set the number travelling
         * to weight. If the edge does exist, add weight to the existing weight if
         * sum_weights is true, otherwise, do nothing if it's false.
         *
         * Return true if the edge was created.
         */

        Edge e;
        bool edge_added;
        std::tie(e, edge_added) = graph.add_edge(vsrc, vdst);

        if (edge_added) {
            weights[e] = weight;
            return;
        } else if (sum_weights) {
            weights.at(e) += weight;
            return;
        } else {
            return;
        }
    }

    auto remove_all_edges() -> void {
        /*
         * Remove all couplings and edges from the network.
         */
        graph.remove_all_edges();
        weights.clear();
    }
public:

};

#endif //EPIGRAPH_META_POP_NETWORK_H
