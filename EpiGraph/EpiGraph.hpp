//
// Created by roryh on 27/10/2020.
//

#ifndef EPIGRAPH_CPP_EPIGRAPH_HPP
#define EPIGRAPH_CPP_EPIGRAPH_HPP

#endif //EPIGRAPH_CPP_EPIGRAPH_HPP
#include <EpiGraph/Core/EpiComp.hpp>
#include <EpiGraph/Core/IO.hpp>
#include <EpiGraph/Core/NetEpiComp.hpp>
#include <EpiGraph/Core/Util.hpp>

#include <EpiGraph/Eigen/EigenUtil.hpp>

#include <EpiGraph/Models/SI.hpp>
#include <EpiGraph/Models/SIR.hpp>
#include <EpiGraph/Models/SIXRD.hpp>

#include <EpiGraph/Random/Distributions.hpp>
#include <EpiGraph/Random/RandomMatrix.hpp>

#include <EpiGraph/Spatial/SpatialUtil.hpp>
#include <EpiGraph/Spatial/TravelStrategies.hpp>

namespace EpiGraph {
    auto print_banner() -> void;

    auto print_banner2() -> void;
}

