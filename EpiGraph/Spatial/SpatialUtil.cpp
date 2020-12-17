//
// Created by roryh on 18/10/2020.
//

#include <EpiGraph/Spatial/SpatialUtil.hpp>
#include <cmath>

namespace EpiGraph {
	auto haversine(double theta) -> double {
		//return (1.0 - cos(theta)) / 2.0;
		return sin(theta / 2.0) * sin(theta / 2.0);
	}

	auto long_lat_distance(double lon1, double lat1, double lon2, double lat2) -> double {
		double pi = 3.14159265;
		return 2 * 6371 * asin(sqrt(haversine((lat2 * pi / 180.0) - (lat1 * pi / 180.0)) +
									cos(lat1 * pi / 180.0) * cos(lat2 * pi / 180.0) *
									haversine((lon2 * pi / 180.0) - (lon1 * pi / 180.0))));
	}

	auto long_lat_distance_2(double lon1, double lat1, double lon2, double lat2) -> double {
		double pi = 3.14;
		return 2 * 6371e3 * asin(sqrt(
				haversine(lon2 - lon1) + cos(lon1) * cos(lon2) * haversine((lat2) - (lat1))));
	}
}