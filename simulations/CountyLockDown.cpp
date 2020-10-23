#include <epigraph/RandomMatrix.hpp>
#include <epigraph/Models/SixrdMetaPop.hpp>
#include <libs/toml11/toml.hpp>

#include <iostream>
#include <random>
#include <chrono>

using namespace Eigen;
using ArrayXXb = Array<bool, Dynamic, Dynamic>;

using Model = SixrdMetaPop;

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

    int i;
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

    // holds the SIXRD state_impl of each node
    Model x(num_nodes);
    x.set_compartment(Model::Sidx, (pop.array() > 0).select(pop.array(), 1));

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
    std::map<std::string, double> county_I = accumulate_groups(x.state().col(Iidx), county);
    std::map<std::string, double> county_N = accumulate_groups(pop, county);


    // holds all the parameters
    const auto &phase_order = toml::find<std::vector<std::string>>(config, "parameters", "order");
    auto current_phase = phase_order[0];
    const auto &phase_params = toml::find(config, "parameters", current_phase);


    x.set_params(Model::beta_idx, VectorXd::Ones(num_nodes) * toml::find<double>(phase_params, "beta"));
    x.set_params(Model::c_idx, VectorXd::Ones(num_nodes) * toml::find<double>(phase_params, "c"));
    x.set_params(Model::mu_idx, VectorXd::Ones(num_nodes) * toml::find<double>(phase_params, "mu"));
    x.set_params(Model::alpha_idx, VectorXd::Ones(num_nodes) * toml::find<double>(phase_params, "alpha"));
    x.set_params(Model::kappa_idx, VectorXd::Ones(num_nodes) * toml::find<double>(phase_params, "kappa"));

    VectorXd max_dist_vec = VectorXd::Ones(num_nodes, 1) * toml::find<double>(phase_params, "max_dist");
    VectorXd compliance_vec = VectorXd::Ones(num_nodes, 1) * toml::find<double>(phase_params, "compliance");
    VectorXi duration_vec = VectorXi::Ones(num_nodes, 1) * toml::find<int>(phase_params, "duration");

    std::vector<int> phase_vec(num_nodes, 0);

    // probability distribution for generating commuter numbers
    ProbDist commuter_dist = ProbDist_from_csv(commuter_distribution_path);

    // Write initial conditions
    write_state(x, "0", full_output_path);
    write_state_totals(x, agg_output_path, false);

    double run_time = 0;
    int max_t = toml::find<int>(config, "parameters", "max_t");
    int max_I = toml::find<int>(config, "parameters", "max_I");

    // responsible for adding edges/travellers to the network
    RandomMatrixGenerator rnd_travel(num_nodes);
    MatrixXd max_dist_mat = general_outer_product(max_dist_vec, max_dist_vec,
                                                  [](double a, double b) { return std::min(a, b); });
    ArrayXXb bool_mat = distance_mat.array() < max_dist_mat.array();
    MatrixXd new_travel_weights = (bool_mat).select(travel_weights.array(),
                                                    travel_weights.array() *
                                                    (1 - compliance_vec.rowwise().replicate(
                                                            travel_weights.cols()).array()));
    rnd_travel.set_row_distributions(new_travel_weights);

    for (int t = 0; t < max_t; t++) {

        std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

        // Compute number of travellers from each node
        VectorXd travel_pop = pop.unaryExpr([&](double x) { return x * commuter_dist(global_engine()); });


        // holds the adjacency matrix
        // Generate movements and store in adj matrix
        Eigen::SparseMatrix<double> adj(num_nodes, num_nodes);
        rnd_travel.distribute_vec_over_matrix_rows(adj, travel_pop);
        x.set_coupling(adj);

        // Update the state_impl matrix
        x.set_state(x.state() + derivative(x));

        // Output to file
        write_state(x, std::to_string(t), full_output_path);
        write_state_totals(x, agg_output_path, true);

        county_I = accumulate_groups(x.state().col(Iidx), county);


        for (int node = 0; node < num_nodes; node++) {
            //decrement duration
            duration_vec[node] -= 1;

            if ((county_I[county[node]] > max_I && phase_vec[node] > 2) ||
                (duration_vec[node] == 0 && county_I[county[node]] > max_I &&
                 phase_vec[node] >= 2)) { // if county gone into lockdown
                phase_vec[node] = 2;
            } else if (duration_vec[node] == 0 &&
                       (phase_vec[node] < phase_order.size() - 1)) { // if node reaches end of phase
                phase_vec[node] += 1;
            } else {
                continue;
            }

            current_phase = phase_order[phase_vec[node]];
            const auto &phase_params = toml::find(config, "parameters", current_phase);

            x.set_params(Model::beta_idx, node, toml::find<double>(phase_params, "beta"));
            x.set_params(Model::c_idx, node, toml::find<double>(phase_params, "c"));
            x.set_params(Model::mu_idx, node, toml::find<double>(phase_params, "mu"));
            x.set_params(Model::alpha_idx, node, toml::find<double>(phase_params, "alpha"));
            x.set_params(Model::kappa_idx, node, toml::find<double>(phase_params, "kappa"));

            max_dist_vec[node] = toml::find<double>(phase_params, "max_dist");
            compliance_vec[node] = toml::find<double>(phase_params, "compliance");
            duration_vec[node] = toml::find<int>(phase_params, "duration");
        }

        max_dist_mat = general_outer_product(max_dist_vec, max_dist_vec,
                                             [](double a, double b) { return std::min(a, b); });
        bool_mat = distance_mat.array() < max_dist_mat.array();
        new_travel_weights = (bool_mat).select(travel_weights.array(),
                                               travel_weights.array() *
                                               (1 - compliance_vec.rowwise().replicate(
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

