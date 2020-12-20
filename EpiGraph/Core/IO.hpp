//
// Created by roryh on 24/06/2020.
//

#ifndef EPIGRAPH_IO_H
#define EPIGRAPH_IO_H

#include <EpiGraph/Core/NetEpiCompBase.hpp>

#include <fstream>
#include <sstream>
#include <vector>

namespace EpiGraph {
template <typename Derived>
auto write_state(const NetEpiCompBase<Derived> &model, const std::string &path,
                 const std::string &header = "") -> void {

  auto &x = model.state();
  std::ofstream myfile;

  myfile.open(path);
  if (!header.empty())
    myfile << header << "\n";

  for (int i = 0; i < x.rows(); i++) {
    for (int j = 0; j < x.cols(); j++) {
      myfile << x(i, j) << ",";
    }
    myfile << x.row(i).sum() << "\n";
  }
  myfile.close();
}
} // namespace EpiGraph

#endif // EPIGRAPH_IO_H
