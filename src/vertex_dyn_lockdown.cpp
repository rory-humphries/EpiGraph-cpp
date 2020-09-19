//
// Created by roryh on 30/06/2020.
//

#include "../include/travel_model.h"
#include "../include/sixrd_model.h"
#include "../include/toml11/toml.hpp"

#include <iostream>
#include <random>
#include <chrono>

struct VData : public VSpatialMetaPop, public VSIXRD, public VTravel {
    VData() {}
};

using EData = EMetaPop;

using MetaPopNetwork = Graph;
using EpiModel = SIXRDModel;
using TravelModel = RandomTravelModel;
using state_type = SIXRDModel::state_type;

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

    // the network structure that holds population and location
    MetaPopNetwork net;

    // maps to hold the network meta data
    std::vector<double> v_pop_map;
    std::vector<SIXRDParam> v_param_map;
    std::vector<std::pair<double, double>> v_pos_map;
    std::unordered_map<Edge, double, EdgeHash> e_pop_map;

    io::CSVReader<3> in(vertex_data);
    in.read_header(io::ignore_extra_column, "long", "lat", "population");

    double lon; double lat; double population;
    while(in.read_row(lon, lat, population)){
        auto v = net.add_vertex();
        v_pop_map.push_back(population);
        v_pos_map.emplace_back(lon, lat);
        v_param_map.emplace_back();
    }

    // responsible for adding edges/travellers to the network
    TravelModel rnd_travel(net.num_vertices());
    rnd_travel.set_commuter_dist(ProbDist_from_csv(commuter_distribution));
    rnd_travel.set_travel_probs(distance_distribution);

    // responsible for the viral dynamics on the network
    EpiModel epi_model(net.num_vertices());
    for (int i = 0; i < epi_model.get_dim(); i++)
        epi_model.set_state(i, 0, v_pop_map[i]);

    // the initially infected vertices
    auto rand_verts = toml::find<std::vector<int>>(config, "parameters", "initial_seed");
    for (auto& rand_v: rand_verts)
        epi_model.infect_vertex(rand_v, 2);

    if (full_output)
        epi_model.write_state("0", full_output_path);
    if (agg_output)
        epi_model.write_compartment_totals(agg_output_path, false);


    int t = 0;
    const auto &phase_order = toml::find<std::vector<std::string>>(config, "parameters", "order");
    size_t phase_idx = 0;
    while (phase_idx < phase_order.size()) {

        const auto &current_phase = phase_order[phase_idx];
        const auto &phase_params = toml::find(config, "parameters", current_phase);

        SIXRDParam sixrd_param{};

        sixrd_param.beta = toml::find<double>(phase_params, "beta");
        sixrd_param.mu = toml::find<double>(phase_params, "mu");
        sixrd_param.c = toml::find<double>(phase_params, "c");
        sixrd_param.alpha = toml::find<double>(phase_params, "alpha");
        sixrd_param.kappa = toml::find<double>(phase_params, "kappa");

        rnd_travel.set_max_distance(toml::find<double>(phase_params, "max_dist"));
        rnd_travel.set_compliance(toml::find<double>(phase_params, "compliance"));

        int duration = toml::find<int>(phase_params, "duration");

        for (int tau = 0; tau < duration; tau++, t++) {
            std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

            // Generate movements
            rnd_travel.add_out_travels(net, v_pos_map, v_pop_map, e_pop_map);

            // Update state
            epi_model.update(net, e_pop_map, sixrd_param);

            // output
            std:: cout << t<< ", " << current_phase << std::endl;

            if (full_output)
                epi_model.write_state(std::to_string(t), full_output_path);
            if (agg_output)
                epi_model.write_compartment_totals(agg_output_path, true);

            epi_model.print_compartment_totals();



            std::cout << net.num_edges() << std::endl;

            std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
            std::cout << time_span.count() << " seconds\n" << std::endl;

            // reset edges and edge maps
            net.remove_all_edges();
            e_pop_map.clear();

            if (epi_model.compertment_total(1) > 1000 && phase_idx > 2) {
                phase_idx = 1;
                break;
            }
        }

        phase_idx++;
    }

    std::cout << "Finished!" << std::endl;

}

