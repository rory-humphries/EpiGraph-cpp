#include "../include/travel_model.h"
#include "../include/sixrd_model.h"
#include "../include/toml11/toml.hpp"

#include <iostream>
#include <random>
#include <chrono>

using MetaPopNetwork = Network<VSpatialMetaPop, EMetaPop>;
using EpiModel = SIXRDModel<Network, VSpatialMetaPop, EMetaPop>;
using TravelModel = RandomTravelModel<MetaPopNetwork>;
using state_type = SIXRDModel<Network, VSpatialMetaPop, EMetaPop>::state_type;

int main() {

    const auto config = toml::parse("../config.toml");

    // input file paths
    const auto &file_paths = toml::find(config, "file_paths");
    const std::string distance_distribution = toml::find<std::string>(file_paths, "travel_distribution");
    const std::string commuter_distribution = toml::find<std::string>(file_paths, "commuter_distribution");
    const std::string vertex_data = toml::find<std::string>(file_paths, "vertex_data");


    // output file paths
    bool agg_output = toml::find<bool>(config, "output", "aggregate");
    std::string agg_output_path;
    if (agg_output)
        agg_output_path = toml::find<std::string>(config, "output", "aggregate_path");

    bool full_output = toml::find<bool>(config, "output", "full");
    std::string full_output_path;
    if (full_output)
        full_output_path = toml::find<std::string>(config, "output", "full_path");


    std::map<int, std::map<std::string, double>> params;
    const auto &phase_order = toml::find<std::vector<std::string>>(config, "parameters", "order");
    int current_phase = 0;
    for (auto &k: phase_order) {
        const auto &phase_params = toml::find(config, "parameters", phase_order[current_phase]);
        params[current_phase]["beta"] = toml::find<double>(phase_params, "beta");
        params[current_phase]["mu"] = toml::find<double>(phase_params, "mu");
        params[current_phase]["c"] = toml::find<double>(phase_params, "c");
        params[current_phase]["alpha"] = toml::find<double>(phase_params, "alpha");
        params[current_phase]["kappa"] = toml::find<double>(phase_params, "kappa");
        params[current_phase]["max_dist"] = toml::find<double>(phase_params, "max_dist");
        params[current_phase]["compliance"] = toml::find<double>(phase_params, "compliance");
        params[current_phase]["duration"] = toml::find<int>(phase_params, "duration");
        current_phase++;
    }


    // the network structure that holds population and location
    MetaPopNetwork net;
    add_metapopulations_from_csv(net, vertex_data);

    // responsible for adding edges/travellers to the network
    TravelModel rnd_travel;
    rnd_travel.set_network(net);
    rnd_travel.travel_dist = ProbDist_from_csv(distance_distribution);
    rnd_travel.commuter_dist = ProbDist_from_csv(commuter_distribution);


    // stores the probability distributions of a traveller going to a destination vertex given the traveller is
    // leaving from vertex v_src. The index of the vector correspond to the vertices of the same index. i.e. vec[v_src]
    std::vector<std::discrete_distribution<Vertex>> travel_probs;
    travel_probs.reserve(net.num_vertices());
    for (Vertex v = 0; v < net.num_vertices(); v++) {
        auto tmp_vec = rnd_travel.travel_probabilities(v, false);
        travel_probs.emplace_back(tmp_vec.begin(), tmp_vec.end());
    }

    // responsible for the viral dynamics on the network
    EpiModel epi_model;
    epi_model.set_network(net);

    // will hold the number belonging to each of S,I,.. for each vertex in the network
    state_type state = epi_model.get_fully_susceptible_state();

    // the initially infected vertices
    auto rand_verts = toml::find<std::vector<int>>(config, "parameters", "initial_seed");

    //auto rand_v = 257;//uni_dist(global_engine());
    for (auto& rand_v: rand_verts)
        epi_model.infect_vertex(state, rand_v, 2);


    if (full_output)
        epi_model.write_state(state, "0", full_output_path);
    if (agg_output)
        epi_model.write_compartment_totals(state, agg_output_path, false);

    int t = 0;
    for (current_phase = 0; current_phase < phase_order.size(); current_phase++) {

        double beta = params.at(current_phase).at("beta");
        double mu = params.at(current_phase).at("mu");
        double c = params.at(current_phase).at("c");
        double alpha = params.at(current_phase).at("alpha");
        double kappa = params.at(current_phase).at("kappa");
        double max_dist = params.at(current_phase).at("max_dist");
        double compliance = params.at(current_phase).at("compliance");
        double duration = params.at(current_phase).at("duration");

        for (int tau = 0; tau < duration; tau++, t++) {


            epi_model.m_beta = beta;
            epi_model.m_mu = mu;
            epi_model.m_alpha = alpha;
            epi_model.m_kappa = kappa;
            epi_model.m_c = c;

            rnd_travel.set_max_dist(max_dist);
            rnd_travel.set_compliance(compliance);


            std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();


            // Generate movements
            for (Vertex v = 0; v < net.num_vertices(); v++) {
                rnd_travel.add_out_travels(v, travel_probs[v]);
                //rnd_travel.add_out_travels(v);
            }

            // find change in state
            state_type dxdt = epi_model.get_zero_state();
            epi_model(state, dxdt, t);

            double isum = 0;
            int count = 0;
            for (auto &k: dxdt) {
                state[count][0] += dxdt[count][0];
                state[count][1] += dxdt[count][1];
                state[count][2] += dxdt[count][2];
                state[count][3] += dxdt[count][3];
                state[count][4] += dxdt[count][4];

                isum += state[count][1];

                count++;
            }

            epi_model.print_compartment_totals(state);
            std::cout << net.num_edges() << std::endl;

            // Reset the network_ptr, send everyone home and update property maps
            net.remove_all_edges();


            // output results to csv
            if (full_output)
                epi_model.write_state(state, std::to_string(t), full_output_path);
            if (agg_output)
                epi_model.write_compartment_totals(state, agg_output_path, true);


            std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(
                    t2 - t1);
            std::cout << time_span.count() << " seconds\n" << std::endl;
        }

    }

    std::cout << "Finished!" <<
              std::endl;

}

