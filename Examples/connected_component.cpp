
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <EpiGraph/Graph/BFS.hpp>
#include <EpiGraph/Graph/Degree.hpp>
#include <fstream>
#include <iostream>
#include <chrono>

using namespace Eigen;
using namespace EpiGraph;

int main(int argc, char *argv[]) {

  std::vector<Triplet<int>> tripletList;

  std::ifstream infile;
  infile.open("../../data/graphs/Slashdot0811.txt");

  std::string line;
  std::string token;
  std::istringstream ss;

  int dim = 0;
  while (getline(infile, line, '\n')) {
    ss = std::istringstream(line);
    std::getline(ss, token, '\t');
    int src = stoi(token);
    std::getline(ss, token, '\t');
    int dst = stoi(token);

    tripletList.emplace_back(dst, src, 1);
    // tripletList.emplace_back(src, dst, 1);

    dim = dim < src ? src : dim;
    dim = dim < dst ? dst : dim;
  }
  infile.close();
  dim++;
  Eigen::SparseMatrix<int> twit_mat(dim, dim);
  twit_mat.setFromTriplets(tripletList.begin(), tripletList.end());
  Eigen::SparseMatrix<int> twit_mat_tr = twit_mat.transpose();

  std::cout << "Edges: " << twit_mat.nonZeros() << std::endl;
  std::cout << "Vertices: " << twit_mat.cols() << std::endl;

  //VectorXi degs = out_deg(twit_mat);
  //std::cout << degs << std::endl;

  std::vector<int> m(100, 0);
  std::chrono::high_resolution_clock::time_point t1 =
          std::chrono::high_resolution_clock::now();
          
  #pragma omp parallel for
  for (int node = 0; node < 100; node += 1) {
    VectorXi v = BFS_parallel(twit_mat, node);
    int n = (v.array()>0).cast<int>().sum();
    m[node] = n;
  }
  std::chrono::high_resolution_clock::time_point t2 =
          std::chrono::high_resolution_clock::now();
      auto time_span =
          std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
  std::cout << time_span.count() << "\n";
  for (auto &k : m)
    std::cout << k << "\n";
}
