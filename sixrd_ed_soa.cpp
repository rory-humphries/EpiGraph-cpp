#include "travel_model.h"
#include "sixrd_model.h"

#include <iostream>
#include <random>
#include <chrono>

double R_0 = 3.7; // reproduction number
double mu = 1.0 / 10.0; // recovery rate
double beta_org = R_0 * mu; // rate of contact and probability of transmission on contact
double beta = R_0 * mu; // rate of contact and probability of transmission on contact
double c = 1.0; // rate of mixing
double alpha = 0.0028; // death rate
double kappa = 0.1; // rate of moving from infected into quarantine



int main() {

    // deals with adding random edges to the netowrk
    RandomTravelModel rnd_travel("../data/processed/ed_soa_travel_prob_mat.csv",
                                 "../data/processed/ed_commuter_probs.csv");


    // The network model
    SIXRDModel network;
    add_metapopulations_from_csv(network, "../data/processed/ed_soa_vertices.csv");

    std::vector<double> tmp_vec(5, 0);
    std::vector<std::vector<double>> state(network.get_graph().num_vertices(), tmp_vec);

    auto rand_v  = 257;//uni_dist(global_engine());
    for (int i = 0; i < network.get_graph().num_vertices(); i++) {
        if (i != rand_v) {
            state[i][0] = network.populations[i];
            state[i][1] = 0;
            state[i][2] = 0;
            state[i][3] = 0;
            state[i][4] = 0;
        }
        else {
            state[i][0] = network.populations[i] - 2;
            state[i][1] = 2;
            state[i][2] = 0;
            state[i][3] = 0;
            state[i][4] = 0;
        }

    }

    std::cout << rand_v << " is the starting vertex" << std::endl;

    network.write_state(state, "0", "output");

    for (int t = 1; t < 600; t++) {
        std::cout << "Day " << t << std::endl;
        double max_dist = 1000;
        double compliance = 1.0;

        int march12 = 50;
        int may18 = march12 + 67;
        int june8 = may18 + 21;
        int june29 = june8 + 21;
        int july20 = june29 + 21;
        int august10 = july20 + 21;


        if (t > march12) { // 12th march, initial lockdown
            max_dist = 2.0;
            compliance = 0.4;
            c = 0.6;
            beta = R_0 * mu * c;
        }
        if (t > march12+10) {
            max_dist = 2.0;
            compliance = 0.75;
            c = 0.4;
            beta = R_0 * mu * c;
        }
        if (t > may18) {
            max_dist = 5.0;
            compliance = 0.75;
            c = 0.3;
            beta = R_0 * mu * c;
        }
        if (t > june8) {
            max_dist = 20.0;
            compliance = 0.75;
            c = 0.4;
            beta = R_0 * mu * c;
        }
        if (t > june29) {
            max_dist = 20.0;
            compliance = 0.75;
            c = 0.5;
            beta = R_0 * mu * c;
        }
        if (t > july20) {
            max_dist = 1000.0;
            compliance = 1;
            c = 0.6;
            beta = R_0 * mu * c;
        }
        if (t > august10) {
            max_dist = 1000;
            compliance = 1;
            c = 0.7;
            beta = R_0 * mu * c;
        }
        if (t > august10+21) {
            max_dist = 1000;
            compliance = 1;
            c = 0.9;
            beta = R_0 * mu * c;
        }

        network.m_beta = beta_org;
        network.m_mu = mu;
        network.m_alpha = alpha;
        network.m_kappa = kappa;
        network.m_c = c;

        rnd_travel.m_max_dist = max_dist;
        rnd_travel.m_compliance = compliance;


        std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();


        // Generate movements
        rnd_travel.add_random_edges(network);

        std::vector<double> tvec(5, 0);
        std::vector<std::vector<double>> dxdt(network.get_graph().num_vertices(), tvec);

        network(state, dxdt, 0);

        int count = 0;
        for (auto &k: dxdt) {
            state[count][0] += dxdt[count][0];
            state[count][1] += dxdt[count][1];
            state[count][2] += dxdt[count][2];
            state[count][3] += dxdt[count][3];
            state[count][4] += dxdt[count][4];
            count++;
        }

        double sumi = 0;
        double sums = 0;
        double sumr = 0;
        double sumx = 0;
        double sumd = 0;
        for (auto &k: state) {
            sums += k[0];
            sumi += k[1];
            sumx += k[2];
            sumr += k[3];
            sumd += k[4];
        }

        std :: cout << "S = " << sums << ", I = " << sumi<< ", R = " << sumr << ", X = " << sumx << ", D = " << sumd << ", N = " << sums+sumi+sumx+sumr+sumd << std::endl;

        // Reset the network, send everyone home and update property maps
        network.remove_all_edges();

        // output results to csv
        network.write_state(state, std::to_string(t), "output");


        std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
        std::cout << time_span.count() << " seconds\n" << std::endl;

    }

    std::cout << "Finished!" << std::endl;

}

