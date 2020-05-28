#include "algos.h"
#include "distributions.h"
#include "spatial_graph.h"
#include "spatial_meta_sirx.h"

#include <iostream>
#include <random>
#include <chrono>

double R_0 = 3.7; // reproduction number
double mu = 1.0 / 10.0; // recovery rate
double beta = R_0 * mu; // rate of contact and probability of transmission on contact
double c = 1.0; // rate of mixing
double alpha = 0.0028; // death rate
double kappa = 0.1; // rate of moving from infected into quarantine

/*
double kappa(double t) {
    // variable rate of moving from infected into quarantine
    // changing compliance with quarantine over time.
    if (t <= 30) {
        return 0.1;
    }
    if (t > 30 && t <= 60) {
        return 0.25;
    }
    else {
        return 0.1;
    }
}
*/
// forward declaration of function to build graph from csv
auto add_edges_from_csv(SpatialGraph& g, std::map<Edge, edata>& emap, std::string path_to_file) -> void;
auto add_vertices_from_csv(SpatialGraph& g, std::map<Vertex, vdata>& vmap, std::string path_to_file) -> void;

int main() {
    std::map<Vertex, vdata> vmap; // holds vertex data
    std::map<Edge, edata> emap; // holds edge data
    SpatialGraph g;
    add_vertices_from_csv(g, vmap, "../data/processed/ed_vertices.csv");
    add_edges_from_csv(g, emap, "../data/processed/ed_edges.csv");

    std::uniform_int_distribution<> uni_dist(0, g.num_vertices());
    auto rand_v  = 2000;
    vmap[rand_v].I = 1/(double)vmap[rand_v].N;
    vmap[rand_v].S = 1 - vmap[rand_v].I;
    rand_v  = 2001;
    vmap[rand_v].I = 1/(double)vmap[rand_v].N;
    vmap[rand_v].S = 1 - vmap[rand_v].I;
    rand_v  = 2002;
    vmap[rand_v].I = 1/(double)vmap[rand_v].N;
    vmap[rand_v].S = 1 - vmap[rand_v].I;

    to_csv(g, vmap, emap, "0", "output");

    for (int t = 1; t < 365; t++) {
        std::cout << "Day " << t << std::endl;
        double max_dist = 1000;
        double compliance = 1.0;

        int march12 = 50;
        int may18 = march12 + 67;
        int june8 = may18 + 21;
        int june29 = june8 + 21;
        int july20 = june29 + 21;
        int august10 = july20 + 21;

        /*
        if (t > march12) { // 12th march, initial lockdown
            max_dist = 2.0;
            compliance = 0.4;
            c = 0.6;
            beta = R_0 * mu * c;
        }
        if (t > march12+10) { // 12th march, initial lockdown
            max_dist = 2.0;
            compliance = 0.8;
            c = 0.4;
            beta = R_0 * mu * c;
        }
        if (t > may18) {
            max_dist = 5.0;
            compliance = 0.8;
            c = 0.3;
            beta = R_0 * mu * c;
        }
        if (t > june8) {
            max_dist = 20.0;
            compliance = 0.8;
            c = 0.5;
            beta = R_0 * mu * c;
        }
        if (t > june29) {
            max_dist = 20.0;
            compliance = 0.8;
            c = 0.6;
            beta = R_0 * mu * c;
        }
        if (t > july20) {
            max_dist = 1000.0;
            compliance = 0;
            c = 0.8;
            beta = R_0 * mu * c;
        }
        if (t > august10) {
            max_dist = 1000;
            compliance = 0;
            c = 1;
            beta = R_0 * mu * c;
        }
         */

        std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

        // Generate movements
        add_edges_commuters(g, vmap, emap, max_dist, compliance);

        // Generate interactions
        gen_interactions(g, vmap, emap, beta, mu, kappa, alpha);

        // Reset the network, send everyone home and update property maps
        reset_movements(g, vmap, emap);
        //gen_interactions(g, vmap, emap, beta, mu, kappa);

        // output results to csv
        to_csv(g, vmap, emap, std::to_string(t), "output");


        print_status(g, vmap);


        std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
        std::cout << time_span.count() << " seconds\n" << std::endl;
    }

    std::cout << "Finished!" << std::endl;

}
auto add_vertices_from_csv(SpatialGraph& g,  std::map<Vertex, vdata>& vmap, std::string path_to_file) -> void {
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

                getline(s, field, ',');
                commuters = stoi(field);
            }
            auto v = g.add_vertex(longatude, latatude);
            std::vector<int> wt_vec(commuters, 0);
            vmap[v].waiting_times = wt_vec;
            vmap[v].N = population;
            vmap[v].S = 1;
            vmap[v].I = 0;
            vmap[v].R = 0;
            vmap[v].X = 0;
        }
    }
    myfile.close();
}
auto add_edges_from_csv(SpatialGraph& g, std::map<Edge, edata>& emap, std::string path_to_file) -> void {
    std::ifstream myfile;
    myfile.open(path_to_file);

    int v_src;
    int v_dst;
    int commuters;
    std::string line;

    if (myfile.is_open()) {
        getline(myfile, line);
        while (getline(myfile, line)) {
            std::istringstream s(line);
            std::string field;



            while (getline(s, field, ',')) {
                v_src = stoi(field);

                getline(s, field, ',');
                v_dst = stoi(field);

                getline(s, field, ',');
                commuters = stoi(field);
            }
            Edge e;
            bool edge_added;
            std::tie(e, edge_added) = g.add_edge(v_src, v_dst);
            emap[e];

            emap[e].commuters = commuters;

            emap[e].N = 0;
            emap[e].S = 1;
            emap[e].I = 0;
            emap[e].R = 0;
            emap[e].X = 0;
        }
    }
    myfile.close();
}
