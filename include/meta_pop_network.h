//
// Created by roryh on 23/06/2020.
//

#ifndef EPIGRAPH_META_POP_NETWORK_H
#define EPIGRAPH_META_POP_NETWORK_H

#include <vector>
#include <map>
#include "edge.h"
#include "graph.h"

struct VMetaPop {
    /*
     * A vertex property that holds a population
     */
    using population_type = double;

    population_type population;
};

struct EMetaPop {
    /*
     * An edge property that holds a population
     */
    using population_type = double;

    population_type population;
};

struct VSpatialMetaPop {
    /*
     * A vertex property that holds a population and location
     */
    using population_type = double;
    using position_type = std::pair<double, double>;

    population_type population;
    position_type position;
};

#endif //EPIGRAPH_META_POP_NETWORK_H
