
#include <EpiGraph/Graph/Degree.hpp>
#include <iostream>

using namespace Eigen;
using namespace EpiGraph;

int main(int argc, char *argv[]) {

  MatrixXi undir_adj = MatrixXi::Zero(4, 4);

  undir_adj(0, 0) = 0;
  undir_adj(0, 1) = 1;
  undir_adj(0, 2) = 1;
  undir_adj(0, 3) = 0;
  undir_adj(1, 1) = 1;
  undir_adj(1, 2) = 1;
  undir_adj(1, 3) = 0;
  undir_adj(2, 2) = 1;
  undir_adj(2, 3) = 1;
  undir_adj(3, 3) = 0;

  std::cout << "undir_adj = \n" << undir_adj << std::endl;

  std::cout << "\nundirected degrees of undir_adj: "
            << deg(undir_adj).transpose() << std::endl;
  std::cout << "\nundirected degree of node 2 of undir_adj: "
            << deg(undir_adj, 2) << std::endl;

  MatrixXi adj = MatrixXi::Zero(3, 3);

  adj(1, 0) = 1;
  adj(2, 0) = 20; // 21

  adj(0, 1) = 1;
  adj(2, 1) = 5; // 6

  adj(0, 2) = 3;
  adj(1, 2) = 0; // 3

  std::cout << "\nadj = \n" << adj << std::endl;

  SparseMatrix<int> adj_sp = adj.sparseView();

  std::cout << "\nadj = \n" << adj_sp << std::endl;

  std::cout << "\nundirected degrees of adj: " << deg(adj).transpose()
            << std::endl;
  std::cout << "\nundirected degrees of adj+adj: " << deg(adj + adj).transpose()
            << std::endl;
  std::cout << "\nUndirected degree of sparse adj: " << deg(adj_sp).transpose()
            << std::endl;
  std::cout << "\nUndirected degree of node 1 of sparse adj: " << deg(adj_sp, 1)
            << std::endl;
  std::cout << "\nOut degrees of adj: " << out_deg(adj).transpose()
            << std::endl;
  std::cout << "\nIn degrees of adj: " << in_deg(adj).transpose() << std::endl;
  std::cout << "\nOut degrees of sparse adj: " << out_deg(adj_sp).transpose()
            << std::endl;
}
