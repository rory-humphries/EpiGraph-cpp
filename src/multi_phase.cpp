#include "../include/travel_model.h"
#include "../include/sixrd_model.h"
#include "../include/toml11/toml.hpp"

#include <iostream>
#include <random>
#include <chrono>


using TravelModel = RandomTravelModelMat;

int main() {

    // toml stuff to read in all configs and paths
    const auto config = toml::parse("../config.toml");

    // input file paths
    const auto &file_paths = toml::find(config, "file_paths");
    const std::string distance_distribution = toml::find<std::string>(file_paths, "travel_distribution");
    const std::string commuter_distribution = toml::find<std::string>(file_paths, "commuter_distribution");
    const std::string vertex_data = toml::find<std::string>(file_paths, "vertex_data");

    // output file paths
    std::string agg_output_path = toml::find<std::string>(config, "output", "aggregate_path");
    std::string full_output_path = toml::find<std::string>(config, "output", "full_path");
    std::string R0_path = toml::find<std::string>(config, "output", "R0_path");

    int dim=0;
    std::ifstream file(vertex_data);
    std::string line;
    while (getline(file, line))
        dim++;
    dim--;

    io::CSVReader<3> in(vertex_data);
    in.read_header(io::ignore_extra_column, "long", "lat", "population");

    // maps to hold the network meta data
    Eigen::Matrix<double, Eigen::Dynamic, 1> pop_vec(dim, 1);
    Eigen::Matrix<double, Eigen::Dynamic, 2> pos_mat(dim, 2);

    double lon, lat, population;
    int i = 0;
    while (in.read_row(lon, lat, population)) {
        pop_vec[i] = population;
        pos_mat(i,0) = lon;
        pos_mat(i,1) = lat;

        i++;
    }


    // responsible for adding edges/travellers to the network
    TravelModel rnd_travel(dim);
    rnd_travel.set_commuter_dist(ProbDist_from_csv(commuter_distribution));
    rnd_travel.set_travel_probs(distance_distribution);

    SIXRD_state<double> x = SIXRD_state<double>::Zero(dim, 5);

    x.col(Sidx) = pop_vec; // set entire population as susceptible

    auto rand_verts = toml::find<std::vector<int>>(config, "parameters", "initial_seed");
    for (auto &rand_v: rand_verts)
        infect_SIXRD_state(x, rand_v, 2);

    write_net_SIXRD_state(x, "0", full_output_path);
    write_net_SIXRD_state_totals(x, agg_output_path, false);


    std::vector<double> R0_vec;

    int t = 0;
    const auto &phase_order = toml::find<std::vector<std::string>>(config, "parameters", "order");
    for (const auto &current_phase : phase_order) {

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
            // output
            std::cout << t << std::endl;

            std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

            // Generate movements and store in adj matrix
            Eigen::SparseMatrix<double> adj(x.rows(), x.rows());
            rnd_travel.add_out_travels(adj, pos_mat, pop_vec);

            // Update the state matrix
            net_SIXRD_ode(x, adj, sixrd_param);

            // Compute the reproduction number
            Eigen::VectorXd x_totals = x.colwise().sum().transpose();
            double R0 = SIXRD_R0(x_totals, sixrd_param);
            R0_vec.push_back(R0);

            // Output to file
            write_net_SIXRD_state(x, std::to_string(t), full_output_path);
            write_net_SIXRD_state_totals(x, agg_output_path, true);

            // Output to console
            std::cout << x.colwise().sum() << std::endl;
            std::cout << "R0 = " << R0 << std::endl;

            std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
            auto time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
            std::cout << time_span.count() << " seconds\n" << std::endl;
        }
        write_vector(R0_vec, R0_path);
    }

    std::cout << "Finished!" << std::endl;

}

