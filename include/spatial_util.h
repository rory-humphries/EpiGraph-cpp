//
// Created by roryh on 22/06/2020.
//

#ifndef EPIGRAPH_SPATIAL_UTIL_H
#define EPIGRAPH_SPATIAL_UTIL_H

#include <utility>
#include <cmath>

auto haversine(double theta) -> double {
    return (1.0 - cos(theta)) / 2.0;
}

auto long_lat_distance(double lon1, double lat1, double lon2, double lat2) -> double {
    double pi = 3.14;
    return 2 * 6371e3 * asin(sqrt(haversine((lon2 * pi / 180.0) - (lon1 * pi / 180.0)) +
                                  cos(lon1 * pi / 180.0) * cos(lon2 * pi / 180.0) *
                                  haversine((lat2 * pi / 180.0) - (lat1 * pi / 180.0))));
}


#endif //EPIGRAPH_SPATIAL_UTIL_H
