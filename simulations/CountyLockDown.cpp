#include <epigraph/Random/RandomMatrix.hpp>
#include <epigraph/Core/NetEpiComp.hpp>
#include <epigraph/epigraph.hpp>
#include <libs/toml11/toml.hpp>

#include <iostream>
#include <random>
#include <chrono>

using namespace Eigen;
using ArrayXXb = Array<bool, Dynamic, Dynamic>;

using Model = NetEpiComp<Eigen::Dynamic>;

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

    int max_t = toml::find<int>(config, "parameters", "max_t");
    int max_I = toml::find<int>(config, "parameters", "max_I");


    // containers to hold the network data

    // holds the population of each node
    auto pop = read_matrix<VectorXd>(node_population_path, true);
    int num_nodes = pop.rows();

    // holds the SIXRD state_impl of each node
    Model x(num_nodes, 5);
    x.set_state(SixrdId::S, (pop.array() > 0).select(pop.array(), 1));

    // Add initial infections
    auto rand_verts = toml::find<std::vector<int>>(config, "parameters", "initial_seed");
    for (auto &rand_v: rand_verts) x.add_infected(rand_v, 2);

    // holds the weights/probabilities of each node travelling to another
    MatrixXd travel_weights = read_matrix<MatrixXd>(travel_weights_path);

    // holds the long lat position of each node
    MatrixX2d pos_mat = read_matrix<MatrixX2d>(node_long_lat_path, true);

    // holds the distance in km between each node
    MatrixXd distance_mat = distance_matrix<MatrixXd>(pos_mat);

    // holds the county of each node
    std::vector<std::string> county = read_vector<std::string>(node_county_path, true);
    std::map<std::string, double> county_I = accumulate_groups(x.state().col(SixrdId::I), county);
    std::map<std::string, double> county_N = accumulate_groups(pop, county);


    // holds all the parameters
    x.set_params(SixrdParamId::beta, VectorXd::Ones(num_nodes) * beta_list[0]);
    x.set_params(SixrdParamId::c, VectorXd::Ones(num_nodes) * c_list[0]);
    x.set_params(SixrdParamId::mu, VectorXd::Ones(num_nodes) * mu_list[0]);
    x.set_params(SixrdParamId::alpha, VectorXd::Ones(num_nodes) * alpha_list[0]);
    x.set_params(SixrdParamId::kappa, VectorXd::Ones(num_nodes) * kappa_list[0]);

    VectorXd cur_max_dist = VectorXd::Ones(num_nodes, 1) * max_dist_list[0];
    VectorXd cur_compliance = VectorXd::Ones(num_nodes, 1) * compliance_list[0];
    VectorXi cur_duration = VectorXi::Ones(num_nodes, 1) * duration_list[0];

    std::vector<int> phase_vec(num_nodes, 0);

    // probability distribution for generating commuter numbers
    ProbDist commuter_dist = ProbDist_from_csv(commuter_distribution_path);

    // Write initial conditions
    write_state(x, "0", full_output_path);
    write_state_totals(x, agg_output_path, false);

    double run_time = 0;


    // responsible for adding edges/travellers to the network
    RandomMatrixGenerator rnd_travel(num_nodes);
    MatrixXd max_dist_mat = general_outer_product(cur_max_dist, cur_max_dist,
                                                  [](double a, double b) { return std::min(a, b); });
    ArrayXXb bool_mat = distance_mat.array() < max_dist_mat.array();
    MatrixXd new_travel_weights = (bool_mat)
            .select(travel_weights.array(),
                    travel_weights.array() * (1 - cur_compliance.rowwise().replicate(travel_weights.cols()).array()));
    rnd_travel.set_row_distributions(new_travel_weights);

    for (int t = 0; t < max_t; t++) {

        std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

        // Compute number of travellers from each node
        VectorXd travel_pop = pop.unaryExpr([&](double x) { return x * commuter_dist(global_engine()); });


        // holds the adjacency matrix
        // Generate movements and store in adj matrix
        x.set_coupling(rnd_travel.distribute_vec_over_matrix_rows(travel_pop));

        // Update the state_impl matrix
        x.set_state(x.state() + sixrd_meta_pop_ode(x));

        // Output to file
        write_state(x, std::to_string(t), full_output_path);
        write_state_totals(x, agg_output_path, true);

        county_I = accumulate_groups(x.state().col(SixrdId::I), county);


        for (int node = 0; node < num_nodes; node++) {
            //decrement duration
            cur_duration[node] -= 1;

            if ((county_I[county[node]] > max_I && phase_vec[node] > 2) ||
                (cur_duration[node] == 0 && county_I[county[node]] > max_I &&
                 phase_vec[node] >= 2)) { // if county gone into lockdown
                phase_vec[node] = 2;
            } else if (cur_duration[node] == 0 &&
                       (phase_vec[node] < phase_order.size() - 1)) { // if node reaches end of phase
                phase_vec[node] += 1;
            } else {
                continue;
            }

            x.set_params(SixrdParamId::beta, node, beta_list[phase_vec[node]]);
            x.set_params(SixrdParamId::c, node, c_list[phase_vec[node]]);
            x.set_params(SixrdParamId::mu, node, mu_list[phase_vec[node]]);
            x.set_params(SixrdParamId::alpha, node, alpha_list[phase_vec[node]]);
            x.set_params(SixrdParamId::kappa, node, kappa_list[phase_vec[node]]);

            cur_max_dist[node] = max_dist_list[phase_vec[node]];
            cur_compliance[node] = compliance_list[phase_vec[node]];
            cur_duration[node] = duration_list[phase_vec[node]];
        }

        max_dist_mat = general_outer_product(cur_max_dist, cur_max_dist,
                                             [](double a, double b) { return std::min(a, b); });
        bool_mat = distance_mat.array() < max_dist_mat.array();
        new_travel_weights = (bool_mat).select(travel_weights.array(),
                                               travel_weights.array() *
                                               (1 - cur_compliance.rowwise().replicate(
                                                       travel_weights.cols()).array()));
        rnd_travel.set_row_distributions(new_travel_weights);

        // Output to console
        std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
        auto time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
        run_time += time_span.count();
        std::cout << "#################################\n";
        std::cout << "Time step : " << t << "\n";
        print_totals(x);
        std::cout << "\nRun time : " << time_span.count() << "s , Total run time : " << run_time << "s\n";
        std::cout << "#################################\n\n";


    }

    std::cout << "Finished!" << std::endl;


}

