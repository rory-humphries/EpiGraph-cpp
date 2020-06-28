//
// Created by roryh on 24/06/2020.
//

#ifndef EPIGRAPH_IO_H
#define EPIGRAPH_IO_H

#include <vector>
#include <fstream>
#include <sstream>

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

#endif //EPIGRAPH_IO_H
