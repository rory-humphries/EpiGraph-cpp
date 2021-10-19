#include <EpiGraph/Core/IO.hpp>
#include <EpiGraph/EigenUtil/Accumulate.hpp>
#include <EpiGraph/EigenUtil/IO.hpp>
#include <EpiGraph/EigenUtil/LinAlg.hpp>
#include <EpiGraph/EpiGraph.hpp>
#include <EpiGraph/Models/CovidSim.hpp>
#include <EpiGraph/Models/SIXRDNetMetaPop.hpp>
#include <EpiGraph/Models/SIXRDNetMetaPop_c.hpp>
#include <EpiGraph/Random/Distributions.hpp>
#include <EpiGraph/Random/RandomMatrix.hpp>
#include <EpiGraph/Spatial/SpatialUtil.hpp>

#include <Eigen/Sparse>

#include <libs/toml11/toml.hpp>

#include <chrono>
#include <iostream>
#include <random>

using namespace Eigen;
using namespace EpiGraph;
using ArrayXXb = Array<bool, Dynamic, Dynamic>;
using Model = SIXRDNetMetaPop_c<MatrixXd, MatrixXd, SparseMatrix<double>>;

int main(int argc, char *argv[]) {

  std::string config_path;

  if (argc < 2) {
    std::cout << "Enter path to config >";
    std::cin >> config_path;
  } else if (argc == 2) {
    config_path = argv[1];
  } else {
    std::cout << "Incorrect number of command line arguments passed. Only the "
                 "config file path is expected.\n";
    return -1;
  }
  print_banner();

  std::cout << "\n\nReading in data..." << std::flush;

  // toml stuff to read in all configs and paths
  const auto config = toml::parse(config_path);

  // input file paths
  const auto &file_paths = toml::find(config, "file_paths");
  const std::string travel_weights_path =
      toml::find<std::string>(file_paths, "travel_distribution");
  const std::string commuter_distribution_path =
      toml::find<std::string>(file_paths, "commuter_distribution");
  const std::string node_long_lat_path =
      toml::find<std::string>(file_paths, "node_long_lat");
  const std::string node_population_path =
      toml::find<std::string>(file_paths, "node_population");
  const std::string commuting_zone_path =
      toml::find<std::string>(file_paths, "commuting_zone");
  const std::string distance_matrix_path =
      toml::find<std::string>(file_paths, "distance_matrix_path");

  // output file paths
  std::string agg_output_path =
      toml::find<std::string>(config, "output", "aggregate_path") +
      "MinImpactLDZones.csv";
  std::string full_output_path =
      toml::find<std::string>(config, "output", "full_path") +
      "MinImpactLDZones/";
  std::string R0_path = toml::find<std::string>(config, "output", "R0_path");

  // parameter lists
  const auto &phase_order =
      toml::find<std::vector<std::string>>(config, "parameters", "order");
  const auto &beta_list =
      toml::find<std::vector<double>>(config, "parameters", "beta_list");
  const auto &c_list =
      toml::find<std::vector<double>>(config, "parameters", "c_list");
  const auto &mu_list =
      toml::find<std::vector<double>>(config, "parameters", "mu_list");
  const auto &alpha_list =
      toml::find<std::vector<double>>(config, "parameters", "alpha_list");
  const auto &kappa_list =
      toml::find<std::vector<double>>(config, "parameters", "kappa_list");
  const auto &compliance_list =
      toml::find<std::vector<double>>(config, "parameters", "compliance_list");
  const auto &max_dist_list =
      toml::find<std::vector<double>>(config, "parameters", "max_dist_list");
  const auto &duration_list =
      toml::find<std::vector<int>>(config, "parameters", "duration_list");

  int max_t = toml::find<int>(config, "parameters", "max_t");
  int max_I = toml::find<int>(config, "parameters", "max_I");

  // containers to hold the network data

  // holds the population of each node
  auto pop = read_matrix<VectorXd>(node_population_path, true);
  int num_nodes = pop.rows();

  // holds the SIXRD state_impl of each node
  Model x(num_nodes);
  x.set_state(SixrdId::S, pop); // check for zero pops

  // Add initial infections
  auto rand_verts =
      toml::find<std::vector<int>>(config, "parameters", "initial_seed");
  for (auto &rand_v : rand_verts)
    x.move_state(SixrdId::S, SixrdId::I, rand_v, 2);

  // holds the weights/probabilities of each node travelling to another
  MatrixXd travel_weights = read_matrix<MatrixXd>(travel_weights_path);

  // holds the long lat position of each node
  // MatrixX2d pos_mat = read_matrix<MatrixX2d>(node_long_lat_path, true);

  // holds the distance in km between each node
  MatrixXd distance_mat = read_matrix<MatrixXd>(distance_matrix_path);

  // holds the county of each node
  std::vector<std::string> zone =
      read_vector<std::string>(commuting_zone_path, true);
  std::map<std::string, double> zone_I =
      accumulate_groups(x.state().col(SixrdId::I), zone);
  std::map<std::string, double> zone_N = accumulate_groups(pop, zone);

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
  write_matrix(x.state(), full_output_path + "0.csv", "S,I,X,R,D,N");
  // write_state_totals(x, agg_output_path, false);

  std::cout << "\nRunning simulation...";
  double run_time = 0;

  // responsible for adding edges/travellers to the network
  RandomAdjMat rnd_travel(num_nodes);
  MatrixXd max_dist_mat =
      general_outer_product(cur_max_dist, cur_max_dist,
                            [](double a, double b) { return std::min(a, b); });
  // ArrayXXb bool_mat = distance_mat.array() < max_dist_mat.array();
  // MatrixXd new_travel_weights = (bool_mat).select(
  //  travel_weights.array(),
  //  travel_weights.array() *
  //    (1 -
  //    cur_compliance.rowwise().replicate(travel_weights.cols()).array()));
  // rnd_travel.set_distributions(new_travel_weights);
  rnd_travel.set_distributions(travel_weights);

  MatrixXd state_history(
      std::accumulate(duration_list.begin(), duration_list.end(), 10), 6);

  VectorXd rowvec = x.state().colwise().sum();

  max_dist_mat =
      general_outer_product(cur_max_dist, cur_max_dist,
                            [](double a, double b) { return std::min(a, b); });

  for (int t = 0; t < max_t; t++) {

    std::chrono::high_resolution_clock::time_point t1 =
        std::chrono::high_resolution_clock::now();

    // Compute number of travellers from each node
    VectorXd travel_pop = pop * 0.6;
    // pop.unaryExpr([&](double x) { return x *
    // commuter_dist(global_engine());
    // });

    // holds the adjacency matrix
    // Generate movements and store in adj matrix
    MatrixXd tmp = rnd_travel.gen_sparse_mat(travel_pop);

    ArrayXXb bool_mat = distance_mat.array() < max_dist_mat.array();
    tmp = (bool_mat).select(
        tmp.array(),
        tmp.array() *
            (1 - cur_compliance.rowwise().replicate(tmp.cols()).array()));
    x.set_coupling(tmp.array().round().matrix().sparseView());

    // Update the state_impl matrix
    MatrixXd new_state = x.state() + x.dXdt();

    ArrayXXb eps_mat =
        100000 * new_state.col(SixrdId::I).array() / pop.array() <
        1.0 / 100000.0;
    VectorXd state_diff = (eps_mat).select(new_state.col(SixrdId::I), 0);

    new_state.col(SixrdId::S) += state_diff;
    new_state.col(SixrdId::I) -= state_diff;

    x.set_state(new_state);

    // std::cout << new_state.row(1000) << std::endl;

    // Output to file
    write_matrix(x.state(), full_output_path + std::to_string(t) + ".csv",
                 "S,I,X,R,D,N");

    zone_I = accumulate_groups(x.state().col(SixrdId::I), zone);

    for (int node = 0; node < num_nodes; node++) {
      // decrement duration
      cur_duration[node] -= 1;

      if ((zone_I[zone[node]] > max_I && phase_vec[node] > 2) ||
          (cur_duration[node] == 0 && zone_I[zone[node]] > max_I &&
           phase_vec[node] >= 2)) { // if county gone into lockdown
        phase_vec[node] = 2;
      } else if (cur_duration[node] == 0 &&
                 (phase_vec[node] <
                  phase_order.size() - 1)) { // if node reaches end of phase
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
    max_dist_mat = general_outer_product(
        cur_max_dist, cur_max_dist,
        [](double a, double b) { return std::min(a, b); });
    /*

        bool_mat = distance_mat.array() < max_dist_mat.array();
        new_travel_weights = (bool_mat).select(
          travel_weights.array(),
          travel_weights.array() *
            (1 -
             cur_compliance.rowwise().replicate(travel_weights.cols()).array()));
        rnd_travel.set_distributions(new_travel_weights);
        */

    // Output to console
    auto comp_vec = x.state().colwise().sum();

    printf("\n\n");
    std::cout << "Time step : " << t;

    std::cout << "\n\nS : " << comp_vec[SixrdId::S];
    std::cout << "\nI : " << comp_vec[SixrdId::I];
    std::cout << "\nX : " << comp_vec[SixrdId::X];
    std::cout << "\nR : " << comp_vec[SixrdId::R];
    std::cout << "\nD : " << comp_vec[SixrdId::D];
    std::cout << "\nTravels : " << tmp.array().round().sum();

    std::chrono::high_resolution_clock::time_point t2 =
        std::chrono::high_resolution_clock::now();
    auto time_span =
        std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
    run_time += time_span.count();

    std::cout << "\n\nTime per loop : " << time_span.count();
    std::cout << "\nTotal run time : " << run_time;
    printf("\n");
    // write_matrix(travels, agg_output_path, "travels");
  }
  // write_matrix(travels, agg_output_path, "travels");

  std::cout << "Finished!" << std::endl;
}
