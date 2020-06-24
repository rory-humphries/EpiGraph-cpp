//
// Created by roryh on 19/06/2020.
//

#include <map>
#include "graph.h"
#include "spatial_utils.h"
#include "meta_pop_network.h"

#ifndef EPIGRAPH_SIRX_NETWORK_H
#define EPIGRAPH_SIRX_NETWORK_H


class SIXRDModel : public SpatialMetaPopNetwork {

public:
    double m_beta; // chance of infection on contact
    double m_c;
    double m_mu; // rate of recovery
    double m_kappa;
    double m_alpha;

    // constructors
    SIXRDModel() : SpatialMetaPopNetwork(),
                   m_beta(1),
                   m_c(1),
                   m_mu(1),
                   m_kappa(1),
                   m_alpha(1) {}

    SIXRDModel(double beta,
               double c,
               double mu,
               double alpha,
               double kappa,
               double compliance,
               double max_dist) : SpatialMetaPopNetwork(),
                                  m_beta(beta),
                                  m_c(c),
                                  m_mu(mu),
                                  m_alpha(alpha),
                                  m_kappa(kappa) {}


    void operator()(const state_type &x, state_type &dxdt, const double t) {
        /*
         * The ode describing the evolution of the network SIXRD system.
         * The paramater x and dxdt should be of the type std::vector<std::vector<double>> where the inner
         * vector is of size 5 such that
         *
         * x[i][0] and dxdt[i][0] corresponds to the variable S
         * x[i][1] and dxdt[i][1] corresponds to the variable I
         * x[i][2] and dxdt[i][2] corresponds to the variable X
         * x[i][3] and dxdt[i][3] corresponds to the variable R
         * x[i][4] and dxdt[i][4] corresponds to the variable D
         */
        for (int v = 0; v < graph.num_vertices(); v++) {

            double new_s = x[v][0];
            double new_i = x[v][1];
            double new_n = populations[v];


            auto it_pair = graph.out_edges(v);
            for (auto it = it_pair.first; it != it_pair.second; it++) {
                double tmp_n = weights.at(*it);
                new_n -= tmp_n;
                new_i -= tmp_n * x[v][1] / populations[v];
                new_s -= tmp_n * x[v][0] / populations[v];
            }
            it_pair = graph.in_edges(v);
            for (auto it = it_pair.first; it != it_pair.second; it++) {
                double tmp_n = weights.at(*it);
                new_n += tmp_n;
                new_i += tmp_n * x[it->src][1] / populations[it->src];
            }

            it_pair = graph.in_edges(v);
            for (auto it = it_pair.first; it != it_pair.second; it++) {
                double tmp_n = weights.at(*it);
                double tmp_s = tmp_n * x[it->src][0] / populations[it->src];

                dxdt[it->src][0] -= m_beta * m_c * tmp_s * new_i / new_n;
                dxdt[it->src][1] += m_beta * m_c * tmp_s * new_i / new_n;

            }

            dxdt[v][0] -= m_beta * m_c * new_s * new_i / new_n; // S
            dxdt[v][1] +=
                    m_beta * m_c * new_s * new_i / new_n - m_mu * x[v][1] - m_alpha * x[v][1] - m_kappa * x[v][1];// I
            dxdt[v][2] += m_kappa * x[v][1] - m_mu * x[v][2] - m_alpha * x[v][2]; // X
            dxdt[v][3] += m_mu * x[v][1] + m_mu * x[v][2];// R
            dxdt[v][4] += m_alpha * x[v][1] + m_alpha * x[v][2];// D

        }
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

        for (int v = 0; v<get_graph().num_vertices(); v++)
        {
            myfile << longlat[v].first << ",";
            myfile << longlat[v].second << ",";
            myfile << x[v][0] <<",";
            myfile << x[v][1] <<",";
            myfile << x[v][2] <<",";
            myfile << x[v][3] <<",";
            myfile << x[v][4] <<",";
            myfile << populations[v] <<"\n";
        }
        myfile.close();
    }
};




auto add_metapopulations_from_csv(SIXRDModel &g, std::string path_to_file) -> void {
    std::ifstream myfile;
    myfile.open(path_to_file);

    int vertex_id;
    double longatude;
    double latatude;
    int population;
    int commuters;
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
            g.add_metapop(population, longatude, latatude);
        }
    }
    myfile.close();
}


#endif //EPIGRAPH_SIRX_NETWORK_H
