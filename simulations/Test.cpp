#include <EpiGraph/Core/IO.hpp>
#include <EpiGraph/EigenUtil/IO.hpp>
#include <EpiGraph/EpiGraph.hpp>
#include <EpiGraph/Models/SIXRDNetMetaPop.hpp>
#include <EpiGraph/Random/Distributions.hpp>
#include <EpiGraph/Random/RandomMatrix.hpp>
#include <EpiGraph/Spatial/SpatialUtil.hpp>

#include <Eigen/Sparse>

#include <toml11/toml.hpp>

#include <chrono>
#include <iostream>

using namespace Eigen;
using namespace EpiGraph;

using Model = SIXRDNetMetaPop<MatrixXd, RowVectorXd, SparseMatrix<double>>;

int main(int argc, char *argv[]) {
  print_banner2();
  Eigen::initParallel();

  int num_nodes = 3;

  // holds the SIXRD state_impl of each node
  Model x(num_nodes);
  x.set_state(SixrdId::S, Vector3d(100, 100, 500));

  // Add initial infections
  x.move_state(SixrdId::S, SixrdId::I, 0, 2);
  x.set_params(SixrdParamId::beta, 0.37);
  x.set_params(SixrdParamId::c, 0.7);
  x.set_params(SixrdParamId::mu, 0.1);
  x.set_params(SixrdParamId::alpha, 0.0028);
  x.set_params(SixrdParamId::kappa, 0.1);

  Eigen::MatrixXd adj(num_nodes, num_nodes);
  // adj(0,1) = 5;
  adj(0, 2) = 20;
  // adj(1,0) = 16;
  // adj(1,2) = 30;
  adj(2, 0) = 50;
  // adj(2,1) = 50;
  x.set_coupling(adj.sparseView());

  std::cout << "\nRunning simulation...";
  double run_time = 0;
  for (int t = 1; t < 100; t++) {
    // Update the state_impl matrix
    x.set_state(x.state() + dXdt(x));

    MatrixXd NGM = next_gen_matrix(x);

    // Output to console
    auto comp_vec = x.state().colwise().sum();
    std::cout << "Time step : " << t;

    std::cout << "\n\n\33[2KS : " << comp_vec[SixrdId::S];
    std::cout << "\n\33[2KI : " << comp_vec[SixrdId::I];
    std::cout << "\n\33[2KX : " << comp_vec[SixrdId::X];
    std::cout << "\n\33[2KR : " << comp_vec[SixrdId::R];
    std::cout << "\n\33[2KD : " << comp_vec[SixrdId::D];

    std::cout << "\n\n\33[2KR0 : " << NGM.eigenvalues().cwiseAbs();

    printf("\n");
  }
  std::cout << "\nFinished!" << std::endl;
}
