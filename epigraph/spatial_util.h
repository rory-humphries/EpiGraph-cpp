//
// Created by roryh on 22/06/2020.
//

#ifndef EPIGRAPH_SPATIAL_UTIL_H
#define EPIGRAPH_SPATIAL_UTIL_H


auto haversine(double theta) -> double;

auto long_lat_distance(double lon1, double lat1, double lon2, double lat2) -> double;

auto long_lat_distance_2(double lon1, double lat1, double lon2, double lat2) -> double;


#endif //EPIGRAPH_SPATIAL_UTIL_H
