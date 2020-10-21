#include <epigraph/random_matrix.hpp>
#include <epigraph/sixrd_net_model_up.hpp>
#include <toml11/toml.hpp>

#include <iostream>
#include <chrono>

using namespace Eigen;
using Model = SIXRDNetModelUP;

int main(int argc, char *argv[]) {
    try {
        std::string config_path;

        if (argc < 2) {
            std::cout << "Enter path to config >";
            std::cin >> config_path;
        } else if (argc == 2) {
            config_path = argv[1];
        } else {
            std::cout << "Incorrect number of command line arguments passed. Only the config file path is expected.\n";
            return -1;
        }
        print_banner2();
        std::cout << std::fixed << std::setprecision(2);

        // toml stuff to read in all configs and paths
        const auto config = toml::parse(config_path);

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
        auto pop = read_matrix<VectorXd>(node_population_path, true);
        int num_nodes = pop.rows();

        // holds the SIXRD state of each node
        Model x(num_nodes);
        x.set_compartment(Model::Sidx, (pop.array() > 0).select(pop.array(), 1));

        // Add initial infections
        auto rand_verts = toml::find<std::vector<int>>(config, "parameters", "initial_seed");
        for (auto &rand_v: rand_verts) x.add_infected(rand_v, 2);

        // holds the weights/probabilities of each node travelling to another
        MatrixXd travel_weights = read_matrix<MatrixXd>(travel_weights_path);
        MatrixXd new_travel_weights = travel_weights;

        // holds the long lat position of each node
        MatrixX2d pos_mat = read_matrix<MatrixX2d>(node_long_lat_path, true);

        // holds the distance in km between each node
        MatrixXd distance_mat = distance_matrix<MatrixXd>(pos_mat);

        // responsible for adding edges/travellers to the network
        RandomMatrixGenerator rnd_travel(num_nodes);

        // probability distribution for generating commuter numbers
        ProbDist commuter_dist = ProbDist_from_csv(commuter_distribution_path);

        // Write initial conditions
        write_state(x, "0", full_output_path);
        write_state_totals(x, agg_output_path, false);

        double run_time = 0;

        int t = 0;
        const auto &phase_order = toml::find<std::vector<std::string>>(config, "parameters", "order");
        for (const auto &current_phase : phase_order) {

            const auto &phase_params = toml::find(config, "parameters", current_phase);

            x.set_params(Model::beta_idx, toml::find<double>(phase_params, "beta"));
            x.set_params(Model::c_idx, toml::find<double>(phase_params, "c"));
            x.set_params(Model::mu_idx, toml::find<double>(phase_params, "mu"));
            x.set_params(Model::alpha_idx, toml::find<double>(phase_params, "alpha"));
            x.set_params(Model::kappa_idx, toml::find<double>(phase_params, "kappa"));
            double max_dist = toml::find<double>(phase_params, "max_dist");
            double compliance = toml::find<double>(phase_params, "compliance");

            new_travel_weights = (distance_mat.array() < max_dist).select(travel_weights,
                                                                          (1 - compliance) * travel_weights);
            rnd_travel.set_row_distributions(new_travel_weights);

            int duration = toml::find<int>(phase_params, "duration");
            for (int tau = 0; tau < duration; tau++, t++) {

                std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

                // Compute number of travellers from each node
                VectorXd travel_pop = pop.unaryExpr(
                        [&](double x) {
                            return x * commuter_dist(global_engine());
                        });

                // holds the adjacency matrix
                // Generate movements and store in adj matrix
                Eigen::SparseMatrix<double> adj(num_nodes, num_nodes);
                rnd_travel.distribute_vec_over_matrix_rows(adj, travel_pop);
                x.set_coupling(adj);

                // Update the state matrix
                x.set_state(x.state() + derivative(x));

                // Output to file
                write_state(x, std::to_string(t), full_output_path);
                write_state_totals(x, agg_output_path, true);

                // Output to console
                std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
                auto time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
                run_time += time_span.count();

                std::cout << "Run time : " << time_span.count() << "s, ";
                print_totals(x);
                std::cout << std::flush;
                std::cout << "\r";
                //std::cout << x.params();
                //std::cout << "R0 : " << R0 << "\n";
                //std::cout << "Run time : " << time_span.count() << "s , Total run time : " << run_time << "s";

            }
        }

        std::cout << "Finished!" << std::endl;
    }
    catch (const std::exception &e) {
        std::cerr << e.what();
        return -1;
    }

}

