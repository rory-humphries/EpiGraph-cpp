//
// Created by roryh on 19/06/2020.
//

#include <map>
#include "graph.h"
#include "spatial_utils.h"
#include "meta_pop_network.h"

#ifndef EPIGRAPH_SIRX_NETWORK_H
#define EPIGRAPH_SIRX_NETWORK_H


template<template<typename, typename> class TNetwork, typename VData, typename EData>
class SIXRDModel  {

public:
    using state_type = std::vector<std::vector<double>>;
    using population_type = typename VData::population_type;

    using NetworkPtr = std::shared_ptr<TNetwork<VData, EData>>;

    double m_beta; // chance of infection on contact
    double m_c;
    double m_mu; // rate of recovery
    double m_kappa;
    double m_alpha;

    TNetwork<VData, EData>* network_ptr;

    // constructors
    SIXRDModel() : m_beta(1),
                   m_c(1),
                   m_mu(1),
                   m_kappa(1),
                   m_alpha(1) {}

    SIXRDModel(double beta,
               double c,
               double mu,
               double alpha,
               double kappa) : m_beta(beta),
                                  m_c(c),
                                  m_mu(mu),
                                  m_alpha(alpha),
                                  m_kappa(kappa) {}

    auto set_network(TNetwork<VData, EData> &g) -> void {
        network_ptr = &g;
    }

    auto operator()(const state_type &x, state_type &dxdt, const double t) -> void{
        /*
         * The ode describing the evolution of the network_ptr SIXRD system.
         * The paramater x and dxdt should be of the type std::vector<std::vector<double>> where the inner
         * vector is of size 5 such that
         *
         * x[i][0] and dxdt[i][0] corresponds to the variable S
         * x[i][1] and dxdt[i][1] corresponds to the variable I
         * x[i][2] and dxdt[i][2] corresponds to the variable X
         * x[i][3] and dxdt[i][3] corresponds to the variable R
         * x[i][4] and dxdt[i][4] corresponds to the variable D
         */
        for (int v = 0; v < (*network_ptr).num_vertices(); v++) {

            double new_s = x[v][0];
            double new_i = x[v][1];
            double new_n = (*network_ptr).vprop[v].population;


            auto it_pair = (*network_ptr).out_edges(v);
            for (auto it = it_pair.first; it != it_pair.second; it++) {
                double tmp_n = (*network_ptr).eprop[*it].population;
                new_n -= tmp_n;
                new_i -= tmp_n * x[v][1] / (*network_ptr).vprop[v].population;
                new_s -= tmp_n * x[v][0] / (*network_ptr).vprop[v].population;
            }
            it_pair = (*network_ptr).in_edges(v);
            for (auto it = it_pair.first; it != it_pair.second; it++) {
                double tmp_n = (*network_ptr).eprop[*it].population;
                new_n += tmp_n;
                new_i += tmp_n * x[it->src][1] / (*network_ptr).vprop[it->src].population;
            }

            it_pair = (*network_ptr).in_edges(v);
            for (auto it = it_pair.first; it != it_pair.second; it++) {
                double tmp_n = (*network_ptr).eprop[*it].population;
                double tmp_s = tmp_n * x[it->src][0] / (*network_ptr).vprop[it->src].population;

                dxdt[it->src][0] -= m_beta * m_c * tmp_s * new_i / new_n;
                dxdt[it->src][1] += m_beta * m_c * tmp_s * new_i / new_n;

            }

            dxdt[v][0] -= m_beta * m_c * new_s * new_i / new_n; // S
            dxdt[v][1] += m_beta * m_c * new_s * new_i / new_n - m_mu * x[v][1] - m_alpha * x[v][1] - m_kappa * x[v][1];// I
            dxdt[v][2] += m_kappa * x[v][1] - m_mu * x[v][2] - m_alpha * x[v][2]; // X
            dxdt[v][3] += m_mu * x[v][1] + m_mu * x[v][2];// R
            dxdt[v][4] += m_alpha * x[v][1] + m_alpha * x[v][2];// D

        }
    }

    auto get_zero_state() -> state_type {
        std::vector<double> tmp(5, 0);
        std::vector<std::vector<double>> state((*network_ptr).num_vertices(), tmp);
        return std::move(state);
    }

    auto get_fully_susceptible_state() -> state_type {
        std::vector<double> tmp(5, 0);
        std::vector<std::vector<double>> state((*network_ptr).num_vertices(), tmp);
        for (int i = 0; i < (*network_ptr).num_vertices(); i++) {

            state[i][0] = (*network_ptr).vprop[i].population;
            state[i][1] = 0;
            state[i][2] = 0;
            state[i][3] = 0;
            state[i][4] = 0;

        }
        return std::move(state);
    }

    auto infect_vertex(state_type &x, Vertex v, population_type N) -> void {
        x[v][0] = (*network_ptr).vprop[v].population - 2;
        x[v][1] = 2;
        x[v][2] = 0;
        x[v][3] = 0;
        x[v][4] = 0;
    }

    auto infect_vertex(state_type &x, Vertex v) -> void {
            x[v][0] = (*network_ptr).vprop[v].population - 2;
            x[v][1] = 2;
            x[v][2] = 0;
            x[v][3] = 0;
            x[v][4] = 0;
    }

    void write_state(const state_type &x, const std::string& id, std::string dir) {
        std::ofstream myfile;
        std::string newd = dir;

        myfile.open ("../"+ dir + "/" + id + ".csv");

        myfile << "x,";
        myfile << "y,";
        myfile << "S,";
        myfile << "I,";
        myfile << "R,";
        myfile << "D,";
        myfile << "X,";
        myfile << "N\n";

        for (int v = 0; v<(*network_ptr).num_vertices(); v++)
        {
            myfile << (*network_ptr).vprop[v].position.first << ",";
            myfile << (*network_ptr).vprop[v].position.second << ",";
            myfile << x[v][0] <<",";
            myfile << x[v][1] <<",";
            myfile << x[v][2] <<",";
            myfile << x[v][3] <<",";
            myfile << x[v][4] <<",";
            myfile << (*network_ptr).vprop[v].population << "\n";
        }
        myfile.close();
    }

    void write_compartment_totals(const state_type &x, const std::string& path_to_file, bool append_file = false) {

        double sumi = 0;
        double sums = 0;
        double sumr = 0;
        double sumx = 0;
        double sumd = 0;
        for (auto &k: x) {
            sums += k[0];
            sumi += k[1];
            sumx += k[2];
            sumr += k[3];
            sumd += k[4];
        }

        std::ofstream myfile;

        if (append_file)
            myfile.open (path_to_file, std::ios_base::app);
        else {
            myfile.open (path_to_file);
            myfile << "S,";
            myfile << "I,";
            myfile << "X,";
            myfile << "R,";
            myfile << "D\n";
        }


        myfile << sums <<",";
        myfile << sumi <<",";
        myfile << sumx <<",";
        myfile << sumr <<",";
        myfile << sumd <<"\n";
        myfile.close();
    }

    auto print_compartment_totals(const state_type& x) {
        double sumi = 0;
        double sums = 0;
        double sumr = 0;
        double sumx = 0;
        double sumd = 0;
        for (auto &k: x) {
            sums += k[0];
            sumi += k[1];
            sumx += k[2];
            sumr += k[3];
            sumd += k[4];
        }

        std :: cout << "S = " << sums << ", I = " << sumi<< ", R = " << sumr << ", X = " << sumx << ", D = " << sumd << ", N = " << sums+sumi+sumx+sumr+sumd << "\n";
    }
};


template <typename TNetwork>
auto add_metapopulations_from_csv(TNetwork &g, std::string path_to_file) -> void {
    std::ifstream myfile;
    myfile.open(path_to_file);

    int vertex_id;
    double longatude;
    double latatude;
    int population;
    std::string line;

    if (myfile.is_open()) {
        getline(myfile, line);
        while (getline(myfile, line)) {
            std::istringstream s(line);
            std::string field;

            while (getline(s, field, ',')) {
                vertex_id = stoi(field);

                getline(s, field, ',');
                longatude = stof(field);

                getline(s, field, ',');
                latatude = stof(field);

                getline(s, field, ',');
                population = stoi(field);
            }
            auto v = g.add_vertex();
            g.vprop[v].population = population;
            g.vprop[v].position = {longatude, latatude};
        }
    }
    myfile.close();
}


#endif //EPIGRAPH_SIRX_NETWORK_H
