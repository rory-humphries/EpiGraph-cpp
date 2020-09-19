//
// Created by roryh on 24/06/2020.
//

#ifndef EPIGRAPH_IO_H
#define EPIGRAPH_IO_H

#include <vector>
#include <fstream>
#include <sstream>
#include <Eigen/Core>
#include "eigen_util.h"

template<typename Derived>
auto write_vector(Eigen::MatrixBase<Derived> vec, std::string fpath) -> void{
    /*
     * Writes a column vector to csv.
     */

    col_vector_assert(vec);

    std::ofstream myfile;
    myfile.open (fpath + ".csv");

    for (int i = 0; i<vec.rows()-1;i++) {
        myfile << vec[i];
        myfile << ",";
    }
    myfile << vec[vec.rows()-1];
}

template<typename T>
auto write_vector(std::vector<T> vec, std::string fpath) -> void{
    /*
     * Writes a column vector to csv.
     */

    std::ofstream myfile;
    myfile.open (fpath + ".csv");

    for (int i = 0; i<vec.size()-1;i++) {
        myfile << vec[i];
        myfile << ",";
    }
    myfile << vec[vec.size()-1];
}

auto matrix_from_csv(std::string path_to_file) -> std::vector<std::vector<double>> {
    std::ifstream myfile;
    myfile.open(path_to_file);

    std::string line;

    std::vector<std::vector<double>> op_vec;

    if (myfile.is_open()) {
        while (getline(myfile, line)) {
            std::istringstream s(line);
            std::string field;

            std::vector<double> line_vec;

            while (getline(s, field, ',')) {
                double x = stof(field);
                line_vec.push_back(x);
            }
            op_vec.push_back(line_vec);
        }
    }
    myfile.close();

    return op_vec;
}

template <typename TNetwork>
void write_edges(TNetwork &x, const std::string& id, std::string dir) {
    std::ofstream myfile;
    std::string newd = dir;

    myfile.open (dir + id + ".csv");

    myfile << "source,";
    myfile << "destination,";
    myfile << "population\n";


    for (int v = 0; v<x.num_vertices(); v++)
    {
        auto it_pair = x.out_edges(v);
        auto it_beg = it_pair.first;
        auto it_end = it_pair.second;

        for (; it_beg!=it_end; it_beg++) {

            myfile << it_beg->src << ",";
            myfile << it_beg->dst << ",";
            myfile << x.eprop[*it_beg].population << "\n";
        }
    }
    myfile.close();
}

#endif //EPIGRAPH_IO_H
