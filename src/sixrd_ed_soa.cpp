#include "../include/travel_model.h"
#include "../include/sixrd_model.h"

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


using MetaPopNetwork = Network<VSpatialMetaPop, EMetaPop>;
using EpiModel = SIXRDModel<Network, VSpatialMetaPop, EMetaPop>;
using TravelModel = RandomTravelModel<MetaPopNetwork>;
using state_type = SIXRDModel<Network, VSpatialMetaPop, EMetaPop>::state_type;

int main() {
    double total_run_time = 0;


    // the network structure that holds population and location
    MetaPopNetwork net;
    add_metapopulations_from_csv(net, "../data/processed/ed_soa_vertices.csv");

    // responsible for adding edges/travellers to the network_ptr
    TravelModel rnd_travel;
    rnd_travel.set_network(net);
    rnd_travel.commuter_dist = ProbDist_from_csv("../data/processed/ed_commuter_probs.csv");
    rnd_travel.travel_dist = ProbDist_from_csv("../data/processed/ed_distance_probs.csv");

    // stores the probability distributions of a traveller going to a destination vertex given the traveller is
    // leaving from vertex v_src. The index of the vector correspond to the vertices of the same index. i.e. vec[v_src]
    std::vector<std::discrete_distribution<Vertex>> travel_probs;
    travel_probs.reserve(net.num_vertices());
    for (Vertex v = 0; v < net.num_vertices(); v++) {
        auto tmp_vec = rnd_travel.travel_probabilities(v, false);
        travel_probs.emplace_back(tmp_vec.begin(), tmp_vec.end());
    }

    // responsible for the viral dynamics on the network
    EpiModel epi_model;
    epi_model.set_network(net);

    // will hold the number belonging to each of S,I,.. for each vertex in the network
    state_type state = epi_model.get_fully_susceptible_state();

    // the initially infected vertex
    auto rand_v = 257;//uni_dist(global_engine());
    epi_model.infect_vertex(state, rand_v);

    // the filename for the output file of the total compartment values
    std::string file_output_path("../" + std::to_string(rand_v) + "2_output.csv");

    epi_model.write_state(state, "0", "output");
    //epi_model.write_compartment_totals(state, file_output_path, false);

    double max_dist = 1000;
    double compliance = 1.0;

    int march12 = 50;
    int may18 = march12 + 67;
    int june8 = may18 + 21;
    int june29 = june8 + 21;
    int july20 = june29 + 21;
    int august10 = july20 + 21;

    double isum = 0;
    int p = 0;
    for (int t = 1; t < 600; t++) {
        std::cout << "Day " << t << std::endl;

        //if (isum > 10000 && t > may18) {
        //    p = may18;
        //}
        if (p > march12) { // 12th march, initial lockdown
            max_dist = 2.0;
            compliance = 0.8;
            c = 0.6;
            beta = R_0 * mu * c;
        }
        if (p > march12 + 10) {
            max_dist = 2.0;
            compliance = 0.8;
            c = 0.4;
            beta = R_0 * mu * c;
        }
        if (p > may18) {
            max_dist = 5.0;
            compliance = 0.8;
            c = 0.3;
            beta = R_0 * mu * c;
        }
        if (p > june8) {
            max_dist = 20.0;
            compliance = 0.8;
            c = 0.4;
            beta = R_0 * mu * c;
        }
        if (p > june29) {
            max_dist = 20.0;
            compliance = 0.8;
            c = 0.5;
            beta = R_0 * mu * c;
        }
        if (p > july20) {
            max_dist = 1000.0;
            compliance = 1;
            c = 0.6;
            beta = R_0 * mu * c;
        }
        if (p > august10) {
            max_dist = 1000;
            compliance = 1;
            c = 0.7;
            beta = R_0 * mu * c;
        }
        if (p > august10 + 21) {
            max_dist = 1000;
            compliance = 1;
            c = 1;
            beta = R_0 * mu * c;
        }
        p++;

        epi_model.m_beta = beta_org;
        epi_model.m_mu = mu;
        epi_model.m_alpha = alpha;
        epi_model.m_kappa = kappa;
        epi_model.m_c = c;

        rnd_travel.set_max_dist(max_dist);
        rnd_travel.set_compliance(compliance);


        std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();


        // Generate movements
        for (Vertex v = 0; v < net.num_vertices(); v++) {
            rnd_travel.add_out_travels(v, travel_probs[v]);
            //rnd_travel.add_out_travels(v);
        }

        state_type dxdt = epi_model.get_zero_state();
        // find rate of change
        epi_model(state, dxdt, t);

        isum = 0;
        int count = 0;
        for (auto &k: dxdt) {
            state[count][0] += dxdt[count][0];
            state[count][1] += dxdt[count][1];
            state[count][2] += dxdt[count][2];
            state[count][3] += dxdt[count][3];
            state[count][4] += dxdt[count][4];

            isum += state[count][1];

            count++;
        }

        epi_model.print_compartment_totals(state);
        std::cout << net.num_edges() << std::endl;

        // Reset the network_ptr, send everyone home and update property maps
        net.remove_all_edges();


        // output results to csv
        epi_model.write_state(state, std::to_string(t), "output");
        //epi_model.write_compartment_totals(state, file_output_path, true);


        std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
        total_run_time += time_span.count();
        std::cout << time_span.count() << " seconds\n" << std::endl;
        std::cout << "Avg loop time = " << total_run_time / t << " seconds\n" << std::endl;

    }

    std::cout << "Finished!" << std::endl;

}

