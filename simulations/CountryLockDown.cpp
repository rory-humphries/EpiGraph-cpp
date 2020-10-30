#include <EpiGraph/EpiGraph.hpp>
#include <toml11/toml.hpp>

#include <iostream>
#include <chrono>

using namespace Eigen;
using namespace EpiGraph;

using Model = NetEpiComp<1>;

int main(int argc, char *argv[]) {
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
        Eigen::initParallel();


        std::cout << "\n\nReading in data..." << std::flush;
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

        // parameter lists
        const auto &phase_order = toml::find<std::vector<std::string>>(config, "parameters", "order");
        const auto &beta_list = toml::find<std::vector<double>>(config, "parameters", "beta_list");
        const auto &c_list = toml::find<std::vector<double>>(config, "parameters", "c_list");
        const auto &mu_list = toml::find<std::vector<double>>(config, "parameters", "mu_list");
        const auto &alpha_list = toml::find<std::vector<double>>(config, "parameters", "alpha_list");
        const auto &kappa_list = toml::find<std::vector<double>>(config, "parameters", "kappa_list");
        const auto &compliance_list = toml::find<std::vector<double>>(config, "parameters", "compliance_list");
        const auto &max_dist_list = toml::find<std::vector<double>>(config, "parameters", "max_dist_list");
        const auto &duration_list = toml::find<std::vector<int>>(config, "parameters", "duration_list");

        // containers to hold the network data

        // holds the population of each node
        auto pop = read_matrix<VectorXd>(node_population_path, true);
        int num_nodes = pop.rows();

        // holds the SIXRD state_impl of each node
        Model x(num_nodes, 5, 5);
        x.set_state(SixrdId::S, (pop.array() > 0).select(pop.array(), 1));

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

        std::cout << "\nRunning simulation...";
        double run_time = 0;
        int t = 0;
        for (int current_phase = 0; current_phase < phase_order.size(); current_phase++) {

            x.set_params(SixrdParamId::beta, beta_list[current_phase]);
            x.set_params(SixrdParamId::c, c_list[current_phase]);
            x.set_params(SixrdParamId::mu, mu_list[current_phase]);
            x.set_params(SixrdParamId::alpha, alpha_list[current_phase]);
            x.set_params(SixrdParamId::kappa, kappa_list[current_phase]);
            double max_dist = max_dist_list[current_phase];
            double compliance = compliance_list[current_phase];
            int duration = duration_list[current_phase];

            new_travel_weights = (distance_mat.array() < max_dist).select(travel_weights,
                                                                          (1 - compliance) * travel_weights);
            rnd_travel.set_row_distributions(new_travel_weights);

            for (int tau = 0; tau < duration; tau++, t++) {
                std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

                // Compute number of travellers from each node
                VectorXd travel_pop = pop.unaryExpr(
                        [&](double x) {
                            return x * commuter_dist(global_engine());
                        });

                // holds the adjacency matrix
                // Generate movements and store in adj matrix
                //Eigen::SparseMatrix<double> adj(num_nodes, num_nodes);
                x.set_coupling(rnd_travel.distribute_vec_over_matrix_rows(travel_pop));

                // Update the state_impl matrix
                x.set_state(x.state() + sixrd_meta_pop_ode(x));

                // Output to file
                write_state(x, std::to_string(t), full_output_path);
                write_state_totals(x, agg_output_path, true);

                //double R0 = net_SIXRD_R0(x.state(), x.coupling(), x.params().transpose());
                //Eigen::SparseMatrix<double> T = sixrd_next_gen_matrix(x.state(), x.coupling(), x.params());
                //std::cout << (double)T.nonZeros() << ", " << (double)(T.rows()*T.cols()) << std::endl;

                // Output to console
                auto comp_vec = x.state().colwise().sum();
                printf("\n\n\33[2K");
                std::cout << "Time step : " << t;

                std::cout << "\n\n\33[2KS : " << comp_vec[SixrdId::S];
                std::cout << "\n\33[2KI : " << comp_vec[SixrdId::I];
                std::cout << "\n\33[2KX : " << comp_vec[SixrdId::X];
                std::cout << "\n\33[2KR : " << comp_vec[SixrdId::R];
                std::cout << "\n\33[2KD : " << comp_vec[SixrdId::D];

                //std::cout << "\n\n\33[2KR0 : " << R0;

                std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
                auto time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
                run_time += time_span.count();

                std::cout << "\n\n\33[2KTime per loop : " << time_span.count();
                std::cout << "\n\33[2KTotal run time : " << run_time;
                printf("\n");
                printf("\x1b[12A");
            }
        }

        std::cout << "\nFinished!" << std::endl;

}

