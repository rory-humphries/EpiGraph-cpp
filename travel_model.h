//
// Created by roryh on 24/06/2020.
//

#ifndef EPIGRAPH_TRAVEL_MODEL_H
#define EPIGRAPH_TRAVEL_MODEL_H

#include <vector>
#include <random>
#include <memory>
#include "distributions.h"
#include "meta_pop_network.h"
#include "io.h"
#include "spatial_utils.h"

template<typename TNetwork>
class RandomTravelModel {
public:

    RandomTravelModel() : commuter_dist(), travel_dist(), m_compliance(0), m_max_dist(0) {}
    RandomTravelModel(std::string commuter_dist_path, std::string travel_dist_path) :
            commuter_dist(ProbDist_from_csv(commuter_dist_path)),
            travel_dist(ProbDist_from_csv(travel_dist_path))  {}

    // setters
    auto set_network(TNetwork& g) -> void {
        network_ptr = &g;
    }
    auto set_compliance(double compliance) -> void {
        m_compliance = compliance;
    }
    auto set_max_dist(double max_dist) -> void {
        m_max_dist = max_dist;
    }

    // getters
    auto get_compliance(double compliance) -> double {
        return m_compliance ;
    }
    auto get_max_dist(double max_dist) -> double {
        return m_max_dist;
    }

    auto travel_probabilities(Vertex v_src, bool normalise = true) -> std::vector<double> {
        /*
         * Return the probabilities that given an individual is travelling they go to vertex v_dst.
         * The index of the returned vector corresponds to the vertex v_dst, i.e. vector[v_dst] is the probability
         * of travelling to v_dst from v_src.
         *
         * Note: The probabilites do not sum to 1,
         */
        std::vector<double> probs;
        probs.reserve((*network_ptr).num_vertices());

        double sum_total = 0;
        for (Vertex v_dst = 0; v_dst < (*network_ptr).num_vertices(); v_dst++) {
            double lon1 = 0, lon2 = 0, lat1 = 0, lat2 = 0;
            std::tie(lon1, lat1) = (*network_ptr).vprop[v_src].position;
            std::tie(lon2, lat2) = (*network_ptr).vprop[v_dst].position;
            double d = long_lat_distance(lon1, lat1, lon2, lat2)/1000.0;
            double tmp = travel_dist.get_prob(d);

            if (d > 0) { // if the distance is 0 then there is no circle
                tmp /= 2*3.14*d; // divide by the circumference of the circle that v_dst lies on (possibly weight by population in the future)
                // TODO: the number of nodes that are a distance d away is not 2*pi*r, will probably have to do some sort of binning of the data.
            }

            sum_total += tmp;
            probs.push_back(tmp);
        }
        if (normalise){
            std::transform(probs.begin(), probs.end(), probs.begin(), [sum_total](double d) -> double {return d/sum_total;});
        }
        return probs;
    }
    template<typename RandomGenerator>
    auto add_out_travels(Vertex v_src, RandomGenerator& travel_distribution) -> void {
        /*
         * Add random edges to the network with a random number of travellers along each edge. The total number of
         * travelers out of v_src is decided by the commuter_dist distribution which gives the proportion of the
         * total population in v_src that will travel.
         *
         * travel distribution is a random distribution that supports the operator() which take a random generator and
         * outputs a destination vertex.
         */

        std::uniform_real_distribution<> uni_dist(0, 1);

        // get the number of people travelling
        double travel_prop = commuter_dist(global_engine());
        int N = (*network_ptr).vprop[v_src].population * travel_prop;


        std::map<Vertex, double> edges_to_add;
        for (int n = 0; n < N; n++) {

            // find the destination vertex from the travel distribution
            Vertex v_dst = travel_distribution(global_engine());
            if (v_dst == v_src)
                continue;

            // if distance is greater than max distance only non compliant will travel
            double lon1 = 0, lon2 = 0, lat1 = 0, lat2 = 0;
            std::tie(lon1, lat1) = (*network_ptr).vprop[v_src].position;
            std::tie(lon2, lat2) = (*network_ptr).vprop[v_dst].position;

            if (long_lat_distance(lon1, lat1, lon2, lat2)/1000.0 >= m_max_dist) {
                double u = uni_dist(global_engine());
                if (u >= m_compliance)
                    continue;
            }
            edges_to_add[v_dst] += 1;
        }
        for (auto &k: edges_to_add) {
            // maybe problems if network_ptr is ever undirected
            auto e_pair = (*network_ptr).add_edge(v_src, k.first);
            Edge e = e_pair.first;

            (*network_ptr).eprop[e].population = k.second;
        }
    }
    auto add_out_travels(Vertex v_src) -> void {
        /*
          * Add random edges to the network with a random number of travellers along each edge. The total number of
          * travellers out of v_src is decided by the commuter_dist distribution which gives the proportion of the
          * total population in v_src that will travel.
          *
          * The probabilities are recomputed on every call. If the network is small enough the probabilities should be
          * stored and passed explicitly.
          */
        auto probs = travel_probabilities(v_src, false);
        std::discrete_distribution<> travel_distribution(probs.begin(), probs.end());

        add_out_travels(v_src, travel_distribution);
    }

private:
    double m_compliance;
    double m_max_dist;

    TNetwork* network_ptr;
public:
    ProbDist commuter_dist;
    ProbDist travel_dist;
};


#endif //EPIGRAPH_TRAVEL_MODEL_H
