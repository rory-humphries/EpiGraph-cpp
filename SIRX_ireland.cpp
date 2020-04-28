#include "algos.h"
#include "distributions.h"
#include "spatial_graph.h"
#include "spatial_meta_sirx.h"

#include <iostream>
#include <random>
#include <chrono>

const double R_0 = 6.1; // reproduction number
const double mu = 1.0 / 14.0; // recovery rate
const double beta = R_0 * mu; // rate of contact and probability of transmission on contact
//const double kappa = 0.01; // rate of moving from infected into quarantine

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

// forward declaration of function to build graph from csv
auto build_from_csv(SpatialGraph& g, std::map<Vertex, vdata>& vmap, std::map<Edge, edata>, std::string path_to_file) -> void;

int main() {
    std::map<Vertex, vdata> vmap; // holds vertex data
    std::map<Edge, edata> emap; // holds edge data
    SpatialGraph g;
    build_from_csv(g, vmap, emap, "../data/pop_data.csv");

    std::cout << "Simulating SIR-X model with parameters\n" << "μ = " << mu << "\nβ = " <<
              beta << "\nκ = " << kappa(0) << "\n\n";

    // distributions
    power_law wt_dist(2.6, 1);
    power_law distance_dist(2.6, 1);
    std::uniform_real_distribution<> uni_dist(0, 1);

    // initialise random seed, adding infected individuals to dublin and cork
    add_infected(g, vmap, 300, 215); // dublin
    add_infected(g, vmap, 300, 215);
    add_infected(g, vmap, 300, 215);
    add_infected(g, vmap, 300, 215);
    add_infected(g, vmap, 300, 215);
    add_infected(g, vmap, 300, 215);
    add_infected(g, vmap, 300, 215);

    add_infected(g, vmap, 150, 70);// cork
    add_infected(g, vmap, 150, 70);
    add_infected(g, vmap, 150, 70);

    print_status(g, vmap);
    std::cout << "\n";

    to_csv(g, vmap, emap, "0", "output");

    for (int t = 1; t < 365; t++) {
        std::cout << "Day " << t << std::endl;

        // initialise waiting times
        update_waiting_times(vmap, wt_dist);

        std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

        // Generate movements
        //gen_movements(g, vmap, emap, distance_dist);
        add_edges_gravity_model(g, vmap, emap, 2.6);

        // Generate interactions
        gen_interactions(g, vmap, emap, beta, mu, kappa(t));

        // Reset the network, send everyone home and update property maps
        reset_movements(g, vmap, emap);

        // output results to csv
        to_csv(g, vmap, emap, std::to_string(t), "output");


        print_status(g, vmap);


        std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
        std::cout << time_span.count() << " seconds\n" << std::endl;
    }

    std::cout << "Finished!" << std::endl;

}

auto build_from_csv(SpatialGraph& g, std::map<Vertex, vdata>& vmap, std::map<Edge, edata>, std::string path_to_file) -> void {
    std::ifstream myfile;
    myfile.open(path_to_file);

    int y_pos;
    int x_pos;
    int population;
    std::string line;

    if (myfile.is_open()) {
        getline(myfile, line);

        while (getline(myfile, line)) {
            std::istringstream s(line);
            std::string field;

            while (getline(s, field, ',')) {
                x_pos = stoi(field);

                getline(s, field, ',');
                y_pos = stoi(field);

                getline(s, field, ',');
                population = stoi(field);
                if (population > 0) {
                    Vertex v = g.add_vertex(x_pos, y_pos);
                    vmap[v];
                    std::vector<int> wt_vec(population, 0);
                    vmap[v].waiting_times = wt_vec;
                    vmap[v].N = population;
                    vmap[v].S = 1;
                    vmap[v].I = 0;
                    vmap[v].R = 0;
                    vmap[v].X = 0;
                }
            }
        }
    }
    myfile.close();
}
