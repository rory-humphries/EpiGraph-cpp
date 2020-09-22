#include "../include/travel_model.h"
#include "../include/sixrd_model.h"
#include "../include/toml11/toml.hpp"

#include <iostream>
#include <random>
#include <chrono>

using namespace Eigen;

//using MatrixX2d = Matrix<double, Dynamic, 2>;

int main() {

    // toml stuff to read in all configs and paths
    const auto config = toml::parse("../config.toml");

    // input file paths
    const auto &file_paths = toml::find(config, "file_paths");
    const std::string travel_weights_path = toml::find<std::string>(file_paths, "travel_distribution");
    const std::string commuter_distribution_path = toml::find<std::string>(file_paths, "commuter_distribution");
    const std::string node_long_lat_path = toml::find<std::string>(file_paths, "node_long_lat");
    const std::string node_population_path = toml::find<std::string>(file_paths, "node_population");
    const std::string node_county_path = toml::find<std::string>(file_paths, "node_county");


    // output file paths
    std::string agg_output_path = toml::find<std::string>(config, "output", "aggregate_path");
    std::string full_output_path = toml::find<std::string>(config, "output", "full_path");
    std::string R0_path = toml::find<std::string>(config, "output", "R0_path");


    // containers to hold the network data

    // holds the population of each node
    VectorXd pop = read_matrix<VectorXd>(node_population_path, true);
    int dim = pop.rows();

    // holds the weights/probabilities of each node travelling to another
    MatrixXd travel_weights = read_matrix<MatrixXd>(travel_weights_path);
    assert(travel_weights.rows() == dim && travel_weights.cols() == dim);

    // holds the long lat position of each node
    MatrixX2d pos_mat = read_matrix<MatrixX2d>(node_long_lat_path, true);
    assert(travel_weights.rows() == dim);

    // holds the distance in km between each node
    MatrixXd distance_mat;
    distance_matrix(distance_mat, pos_mat);

    // holds the SIXRD state of each node
    SIXRD_state<double> x = SIXRD_state<double>::Zero(dim, 5);
    x.col(Sidx) = pop;

    // holds the R0 value for each time step
    std::vector<double> R0_vec;

    // responsible for adding edges/travellers to the network
    RandomMatrix rnd_travel(dim);

    // probability distribution for generating commuter numbers
    ProbDist commuter_dist = ProbDist_from_csv(commuter_distribution_path);

    // Add initial infections
    auto rand_verts = toml::find<std::vector<int>>(config, "parameters", "initial_seed");
    for (auto &rand_v: rand_verts)
        infect_SIXRD_state(x, rand_v, 2);

    // Write initial conditions
    write_net_SIXRD_state(x, "0", full_output_path);
    write_net_SIXRD_state_totals(x, agg_output_path, false);

    double run_time = 0;

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

        double max_dist = toml::find<double>(phase_params, "max_dist");
        double compliance = toml::find<double>(phase_params, "compliance");

        MatrixXd new_travel_weights = (distance_mat.array() < max_dist).select(travel_weights,
                                                                               compliance * travel_weights);
        rnd_travel.set_entry_probs(travel_weights);

        int duration = toml::find<int>(phase_params, "duration");
        for (int tau = 0; tau < duration; tau++, t++) {

            std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

            // Compute number of travellers from each node
            VectorXd travel_pop = pop.unaryExpr(
                    [&](double x) {
                        return x * commuter_dist(global_engine());
                    });

            Eigen::SparseMatrix<double> adj(x.rows(), x.rows());
            // Generate movements and store in adj matrix
            rnd_travel.generate(adj, travel_pop);

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
            std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
            auto time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
            run_time += time_span.count();
            RowVectorXd op_vec = x.colwise().sum();
            std::cout << "#################################\n";
            std::cout << "Time step : " << t << "\n";
            std::cout << "S : " << op_vec[Sidx];
            std::cout << ", I : " << op_vec[Iidx];
            std::cout << ", X : " << op_vec[Xidx];
            std::cout << ", R : " << op_vec[Ridx];
            std::cout << ", D : " << op_vec[Didx] << "\n";
            std::cout << "R0 : " << R0 << "\n";
            std::cout << "Run time : " << time_span.count() << "s , Total run time : " << run_time << "s\n";
            std::cout << "#################################\n\n";


        }
        write_vector(R0_vec, R0_path);
    }

    std::cout << "Finished!" << std::endl;

}
