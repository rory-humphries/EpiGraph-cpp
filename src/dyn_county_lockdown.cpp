#include "../include/travel_model.h"
#include "../include/sixrd_model.h"
#include "../include/toml11/toml.hpp"

#include <iostream>
#include <random>
#include <chrono>

using namespace Eigen;

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

    // holds the county of each node
    std::vector<std::string> county = read_vector<std::string>(node_county_path, true);
    std::map<std::string, double> county_I;
    for (auto &k : county) {
        county_I[k] = 0;
    }

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
    SIXRDState<double> x = SIXRDState<double>::Zero(dim, 5);
    x.col(Sidx) = pop;

    // holds all the parameters
    Matrix<double, Dynamic, 5> params = Matrix<double, Dynamic, 5>::Ones(x.rows(), 5);
    const auto &phase_order = toml::find<std::vector<std::string>>(config, "parameters", "order");
    auto current_phase = phase_order[0];
    const auto &phase_params = toml::find(config, "parameters", current_phase);


    params.col(beta_idx) *= toml::find<double>(phase_params, "beta");
    params.col(mu_idx) *= toml::find<double>(phase_params, "mu");
    params.col(c_idx) *= toml::find<double>(phase_params, "c");
    params.col(alpha_idx) *= toml::find<double>(phase_params, "alpha");
    params.col(kappa_idx) *= toml::find<double>(phase_params, "kappa");

    VectorXd max_dist_vec = VectorXd::Ones(x.rows(), 1) * toml::find<double>(phase_params, "max_dist");
    VectorXd compliance_vec = VectorXd::Ones(x.rows(), 1) * toml::find<double>(phase_params, "compliance");
    VectorXi duration_vec = VectorXi::Ones(x.rows(), 1) * toml::find<int>(phase_params, "duration");

    std::vector<int> phase_vec(dim, 0);

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
    int max_t = toml::find<int>(config, "parameters", "max_t");
    int max_I = toml::find<int>(config, "parameters", "max_I");

    MatrixXd mat1 = max_dist_vec.asDiagonal() * MatrixXd::Ones(dim, dim);
    MatrixXd mat2 = MatrixXd::Ones(dim, dim) * max_dist_vec.asDiagonal();

    MatrixXd new_travel_weights = (distance_mat.array() < mat1.array().min(mat2.array())).select(travel_weights,
                                                                                                 compliance_vec *
                                                                                                 travel_weights);
    rnd_travel.set_entry_probs(travel_weights);

    for (int t = 0; t < max_t; t++) {

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
        net_SIXRD_ode_inhom(x, adj, params);;

        // Output to file
        write_net_SIXRD_state(x, std::to_string(t), full_output_path);
        write_net_SIXRD_state_totals(x, agg_output_path, true);

        for (auto &k : county_I) {
            k.second = 0;
        }
        for (int node = 0; node < x.rows(); node++) {
            county_I[county[node]] += x(node, Iidx);
        }
        for (int node = 0; node < x.rows(); node++) {
            //decrement duration
            duration_vec[node]--;

            if (county_I[county[node]] > max_I && phase_vec[node] > 2) { // if county gone into lockdown
                phase_vec[node] = 2;
            } else if (duration_vec[node] == 0 &&
                       (phase_vec[node] < phase_order.size() - 1)) { // if node reaches end of phase
                phase_vec[node]++;
            } else {
                continue;
            }

            auto current_phase = phase_order[phase_vec[node]];
            const auto &phase_params = toml::find(config, "parameters", current_phase);

            params(node, beta_idx) = toml::find<double>(phase_params, "beta");
            params(node, mu_idx) = toml::find<double>(phase_params, "mu");
            params(node, c_idx) = toml::find<double>(phase_params, "c");
            params(node, alpha_idx) = toml::find<double>(phase_params, "alpha");
            params(node, kappa_idx) = toml::find<double>(phase_params, "kappa");

            max_dist_vec[node] = toml::find<double>(phase_params, "max_dist");
            compliance_vec[node] = toml::find<double>(phase_params, "compliance");
            duration_vec[node] = toml::find<int>(phase_params, "duration");
        }

        MatrixXd mat1 = max_dist_vec.asDiagonal() * MatrixXd::Ones(x.rows(), x.rows());
        MatrixXd mat2 = MatrixXd::Ones(x.rows(), x.rows()) * max_dist_vec.asDiagonal();

        MatrixXd new_travel_weights = (distance_mat.array() < mat1.array().min(mat2.array())).select(travel_weights,
                                                                                                     compliance_vec *
                                                                                                     travel_weights);
        rnd_travel.set_entry_probs(travel_weights);

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
        //std::cout << "R0 : " << R0 << "\n";
        std::cout << "Run time : " << time_span.count() << "s , Total run time : " << run_time << "s\n";
        std::cout << "#################################\n\n";


    }

    std::cout << "Finished!" << std::endl;

}

