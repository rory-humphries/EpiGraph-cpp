//
// Created by roryh on 06/04/2020.
//

#ifndef EPIGRAPH_ALGOS_H
#define EPIGRAPH_ALGOS_H

#include <utility>
#include <set>
#include <vector>
#include <algorithm>
#include <iostream>


template <typename OutputIter>
OutputIter bresenhams_circle(double radius, int x0, int y0, OutputIter result) {
    /*
     * Bresenhams circle algorithm / midpoint algorithm. Outputs an integer raster circle of radius centred at position
     * (x0, y0).
     */

    // The circle algo gives the original cell as a circle of radius 1
    // so take the radius as radius = radius + 1

    int x = 0, y = radius;
    int d = 3 - 2 * radius;
    std::vector<int> x_vals;
    std::vector<int> y_vals;

    x_vals.push_back(x);
    y_vals.push_back(y);
    while (y > x) {
        // for each pixel we will
        // draw all eight pixels

        x++;

        // check for decision parameter
        // and correspondingly
        // update d, x, y
        if (d > 0) {
            y--;
            d = d + 4 * (x - y) + 10;
        } else
            d = d + 4 * x + 6;
        x_vals.push_back(x);
        y_vals.push_back(y);
    }

    std::vector<std::pair<int,int>> circle;
    size_t max_size = x_vals.size();
    // go through reflection of each octant
    for (size_t i = 0; i < max_size; i++) {
        *result = std::make_pair(x_vals[i]+x0, y_vals[i]+y0);;
        result++;
    }
    for (size_t i = max_size - 2; i < max_size; i--) {
        *result =  std::make_pair(y_vals[i]+x0, x_vals[i]+y0);
        result++;
    }
    for (size_t i = 1; i < max_size; i++) {
        *result = std::make_pair(y_vals[i]+x0, -x_vals[i]+y0);
        result++;
    }
    for (size_t i = x_vals.size() - 2; i < x_vals.size(); i--) {
        *result = std::make_pair(x_vals[i]+x0, -y_vals[i]+y0);
        result++;
    }
    for (size_t i = 1; i < max_size; i++) {
        *result = std::make_pair(-x_vals[i]+x0, -y_vals[i]+y0);
        result++;
    }
    for (size_t i = max_size - 2; i < max_size; i--) {
        *result = std::make_pair(-y_vals[i]+x0, -x_vals[i]+y0);
        result++;
    }
    for (size_t i = 1; i < max_size; i++) {
        *result = std::make_pair(-y_vals[i]+x0, x_vals[i]+y0);
        result++;
    }
    for (size_t i = max_size - 2; i < max_size; i--) {
        *result = std::make_pair(-x_vals[i]+x0, y_vals[i]+y0);
        result++;
    }

    return result;
}

#endif //EPIGRAPH_ALGOS_H
