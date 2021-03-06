//
// Created by roryh on 02/11/2020.
//

#ifndef EPIGRAPH_CPP_IO_HPP
#define EPIGRAPH_CPP_IO_HPP

#include <EpiGraph/EigenUtil/StaticAsserts.hpp>

#include <Eigen/Dense>

#include <fstream>
#include <sstream>
#include <vector>

namespace EpiGraph {
template <typename Derived>
auto write_vector(Eigen::MatrixBase<Derived> vec, std::string fpath) -> void {
  /*
   * Writes a column vector to csv.
   */

  col_vector_assert(vec);

  std::ofstream myfile;
  myfile.open(fpath);
  if (myfile.fail()) {
    throw std::runtime_error("Failed to open file " + fpath);
  }

  for (int i = 0; i < vec.rows() - 1; i++) {
    myfile << vec[i];
    myfile << ",";
  }
  myfile << vec[vec.rows() - 1];
}

template <typename T>
auto write_vector(std::vector<T> vec, std::string fpath) -> void {
  /*
   * Writes a column vector to csv.
   */

  std::ofstream myfile;
  myfile.open(fpath);
  if (myfile.fail()) {
    throw std::runtime_error("Failed to open file " + fpath);
  }

  for (int i = 0; i < vec.size() - 1; i++) {
    myfile << vec[i];
    myfile << ",";
  }
  myfile << vec[vec.size() - 1];
}

template <typename T>
auto read_vector(std::string path_to_file, bool header = false)
    -> std::vector<T> {
  std::ifstream infile;
  infile.open(path_to_file);
  if (infile.fail()) {
    throw std::runtime_error("Failed to open file " + path_to_file);
  }

  std::string line;

  // If there is a header ignore it
  if (header)
    getline(infile, line, '\n');
  std::vector<T> op_vec;
  while (getline(infile, line, '\n')) {
    op_vec.push_back(line);
  }
  infile.close();

  return op_vec;
}

template <typename Mat>
auto read_matrix(std::string path_to_file, bool header = false) -> Mat {
  using Scalar = typename Mat::Scalar;

  std::ifstream infile;
  infile.open(path_to_file);
  if (infile.fail()) {
    throw std::runtime_error("Failed to open file " + path_to_file);
  }

  std::string line;
  std::vector<Scalar> op_vec;

  // If there is a header ignore it
  if (header)
    getline(infile, line, '\n');

  // Ensure consistent col length
  getline(infile, line, '\n');
  std::vector<Scalar> tmp;
  std::istringstream ss(line);
  std::string token;

  int cols = 0;
  while (std::getline(ss, token, ',')) {
    op_vec.push_back(stod(token));
    cols++;
  }

  int rows = 1;
  while (getline(infile, line, '\n')) {
    ss = std::istringstream(line);
    int curr_cols = 0;
    while (std::getline(ss, token, ',')) {
      op_vec.push_back(stod(token));
      curr_cols++;
    }
    if (curr_cols != cols)
      throw std::invalid_argument("Inconsistent line length");
    rows++;
  }
  infile.close();

  return Eigen::Map<Mat, 0, Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>>(
      op_vec.data(), rows, cols,
      Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>(1, cols));
}

template <typename Derived>
auto write_matrix(const Eigen::MatrixBase<Derived> &x, const std::string &path,
                  const std::string &header = "") -> void {

  std::ofstream myfile;

  myfile.open(path);
  if (!header.empty())
    myfile << header << "\n";

  for (int i = 0; i < x.rows(); i++) {
    for (int j = 0; j < x.cols(); j++) {
      myfile << x(i, j) << ",";
    }
    myfile << "\n";
  }
  myfile.close();
}
} // namespace EpiGraph

#endif // EPIGRAPH_CPP_IO_HPP
