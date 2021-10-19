#include <Eigen/Sparse>
#include <EpiGraph/Core/IO.hpp>
#include <EpiGraph/EigenUtil/Accumulate.hpp>
#include <EpiGraph/EigenUtil/IO.hpp>
#include <EpiGraph/EigenUtil/LinAlg.hpp>
#include <EpiGraph/EpiGraph.hpp>
#include <EpiGraph/Models/CovidSim.hpp>
#include <EpiGraph/Models/SIXRDNetMetaPop.hpp>
#include <EpiGraph/Random/Distributions.hpp>
#include <EpiGraph/Random/RandomMatrix.hpp>
#include <EpiGraph/Spatial/SpatialUtil.hpp>
#include <chrono>
#include <iostream>
#include <libs/toml11/toml.hpp>
#include <random>

using namespace Eigen;
using namespace EpiGraph;

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

  CovidSimParams params;

  // holds the population of each node
  params.population = read_matrix<VectorXd>(node_population_path, true);
  params.dim = params.population.rows();
  params.phase_order =
      toml::find<std::vector<std::string>>(config, "parameters", "order");
  params.num_phases = params.phase_order.size();

  std::vector<PhaseParams> phase_params_list(params.num_phases);

  // parameter lists

  std::vector<double> beta_list =
      toml::find<std::vector<double>>(config, "parameters", "beta_list");
  std::vector<double> c_list =
      toml::find<std::vector<double>>(config, "parameters", "c_list");
  std::vector<double> mu_list =
      toml::find<std::vector<double>>(config, "parameters", "mu_list");
  std::vector<double> alpha_list =
      toml::find<std::vector<double>>(config, "parameters", "alpha_list");
  std::vector<double> kappa_list =
      toml::find<std::vector<double>>(config, "parameters", "kappa_list");
  std::vector<double> compliance_list =
      toml::find<std::vector<double>>(config, "parameters", "compliance_list");
  std::vector<double> max_dist_list =
      toml::find<std::vector<double>>(config, "parameters", "max_dist_list");
  std::vector<int> duration_list =
      toml::find<std::vector<int>>(config, "parameters", "duration_list");

  for (int i = 0; i < params.num_phases; i++) {
    phase_params_list[i].sixrd_params.beta = beta_list[i];
    phase_params_list[i].sixrd_params.c = c_list[i];
    phase_params_list[i].sixrd_params.mu = mu_list[i];
    phase_params_list[i].sixrd_params.alpha = alpha_list[i];
    phase_params_list[i].sixrd_params.kappa = kappa_list[i];

    phase_params_list[i].compliance = compliance_list[i];
    phase_params_list[i].max_dist = max_dist_list[i];
    phase_params_list[i].duration = duration_list[i];
  }

  // set zones which each node belongs to, can be counties etc..
  std::vector<std::string> zone =
      read_vector<std::string>(commuting_zone_path, true);

  // set travel distributions
  MatrixXd travel_weights = read_matrix<MatrixXd>(travel_weights_path);
  RandomAdjMat rnd_travel(params.dim);
  rnd_travel.set_distributions(travel_weights);
  params.random_adj_mat = rnd_travel;

  // set distance matrix
  params.distance_mat = read_matrix<MatrixXd>(distance_matrix_path);

  params.phase_params_list = phase_params_list;
  CovidSim simulation(params);

  // Add initial infections
  auto rand_verts =
      toml::find<std::vector<int>>(config, "parameters", "initial_seed");
  for (auto &rand_v : rand_verts) simulation.add_infections(rand_v, 2);

  // Write initial conditions
  write_matrix(simulation.get_model().state(), full_output_path + "0.csv",
               "S,I,X,R,D,N");

  std::cout << "\nRunning simulation...";
  double run_time = 0;

  for (int t = 0; t < 600; t++) {
    std::chrono::high_resolution_clock::time_point t1 =
        std::chrono::high_resolution_clock::now();

    simulation.update_travels();
    simulation.update_state();
    simulation.update_phase();

    printf("\n\n");
    std::cout << "Time step : " << t;

    simulation.print_compartment_totals();

    std::chrono::high_resolution_clock::time_point t2 =
        std::chrono::high_resolution_clock::now();
    auto time_span =
        std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
    run_time += time_span.count();

    std::cout << "\n\nTime per loop : " << time_span.count();
    std::cout << "\nTotal run time : " << run_time;
    printf("\n");
  }

  std::cout << "Finished!" << std::endl;
}
