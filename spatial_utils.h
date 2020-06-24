//
// Created by roryh on 22/06/2020.
//

#ifndef EPIGRAPH_SPATIAL_UTILS_H
#define EPIGRAPH_SPATIAL_UTILS_H

#include <utility>
#include <cmath>

auto long_lat_distance(std::pair<double, double> v1, std::pair<double, double> v2) -> double {
    double pi = 3.14;
    double lat1 = v1.second;
    double lon1 = v1.first;
    double lat2 = v2.second;
    double lon2 = v2.first;

    double R = 6371e3; // metres
    double phi1 = lat1 * pi / 180.0; // φ, λ in radians
    double phi2 = lat2 * pi / 180.0;
    double delphi = (lat2 - lat1) * pi / 180.0;
    double dellam = (lon2 - lon1) * pi / 180.0;

    double a = sin(delphi / 2) * sin(delphi / 2) +
               cos(phi1) * cos(phi2) *
               sin(dellam / 2) * sin(dellam / 2);
    double c = 2 * atan2(sqrt(a), sqrt(1 - a));

    double d = R * c; // in metres
    return d/1000;
}

#endif //EPIGRAPH_SPATIAL_UTILS_H
