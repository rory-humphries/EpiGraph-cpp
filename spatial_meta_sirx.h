//
// Created by roryh on 05/04/2020.
//

#ifndef EPIGRAPH_SPATIAL_META_SIR_H
#define EPIGRAPH_SPATIAL_META_SIR_H


#include "distributions.h"
#include "spatial_graph.h"


struct vdata {
    std::vector<int> waiting_times;
    double S; // susceptible
    double I; // infected
    double R; // recovered
    double D; // dead
    double X; // quarantined

    double N;


};
struct edata {
    double S;
    double I;
    double R;
    double D;
    double X;

    double N;
    int commuters;
};

auto add_edges_gravity_model(SpatialGraph &g, std::map<Vertex, vdata>& vmap, std::map<Edge, edata>& emap, double d_alpha) -> void {
    /*
     *  Used rejection sampling by assuming sample vertices are chosen according to a power law distribution dij^{-alpha},
     *  then we can sample from the gravity model (mi/mj)*dij^{-alpha} by accepting a sample if for a sample u~u(0,1)
     *  u < (mj/mi), note the fraction is flipped. mi is the population of vertex i, mj the population of vertex j.
     *
     *  vmap and emap are maps that associate each vertex and edge with SIRX data. d_alpha is the exponent of the
     *  power law from which the random distance is drawn.
     */

    int num_travels = 0;
    std::uniform_real_distribution<> uni_dist(0, 2.6);
    std::uniform_int_distribution<> uni_int_dist(0, g.num_vertices());
    power_law distance_dist(d_alpha, 2.6);

    double M, u, d, ang;
    Vertex v_src, v_dst;
    bool rejected_sample;
    for (v_src = 0; v_src < g.num_vertices(); v_src++) {
        for (int n = 0; n < vmap[v_src].N; n++) {
            if (vmap[v_src].waiting_times[n] != 0) {
                continue;
            }
            rejected_sample = true;
            while(rejected_sample) {
                u = uni_dist(global_engine());
                d = distance_dist(global_engine());
                ang = uni_dist(global_engine());

                double new_x = g.x_pos(v_src) + d * cos(2 * 3.14 * ang);
                double new_y = g.y_pos(v_src) + d * sin(2 * 3.14 * ang);
                v_dst = g.closest_vertex(new_x, new_y);

                //if (pow(new_x - g.x_pos(v_dst), 2) + pow(new_y - g.y_pos(v_dst), 2) > 1 || v_src == v_dst)
                if (pow(new_x - g.x_pos(v_dst), 2) + pow(new_y - g.y_pos(v_dst), 2) > 10 || v_src == v_dst) {
                    continue;
                }
                M = (double)vmap[v_src].N / (double) vmap[v_dst].N;

                if (u < M) {
                    rejected_sample = false;
                    num_travels++;

                    Edge e;
                    bool edge_added;
                    std::tie(e, edge_added) = g.add_edge(v_src, v_dst);

                    if (!edge_added) {
                        emap[e].N += 1;
                        vmap[v_src].N -= 1;
                    } else {
                        emap[e];

                        vmap[v_src].N -= 1;
                        emap[e].N += 1;
                        emap[e].S = vmap[v_src].S;
                        emap[e].I = vmap[v_src].I;
                        emap[e].R = vmap[v_src].R;
                        emap[e].X = vmap[v_src].X;
                    }
                }
            }
        }
    }
    std::cout << "Number of travels = " << num_travels << std::endl;
}

auto add_edges_commuters(SpatialGraph &g, std::map<Vertex, vdata>& vmap, std::map<Edge, edata>& emap, double max_dist, double compliance) -> void {

    Vertex v_src, v_dst;
    for (v_src = 0; v_src < g.num_vertices(); v_src++) {
        auto it_pair = g.out_edges(v_src);
        for(; it_pair.first!=it_pair.second; it_pair.first++) {
            auto e = *it_pair.first; // edge
            if (g.long_lat_distance(v_src, e.dst)/1000 >= max_dist) {
                emap[e].N = (int)(1-compliance)*emap[e].commuters;
                vmap[v_src].N -= emap[e].N;
                emap[e].S = vmap[v_src].S;
                emap[e].I = vmap[v_src].I;
                emap[e].R = vmap[v_src].R;
                emap[e].X = vmap[v_src].X;
            }
            else {
                emap[e].N = emap[e].commuters;
                vmap[v_src].N -= emap[e].commuters;
                emap[e].S = vmap[v_src].S;
                emap[e].I = vmap[v_src].I;
                emap[e].R = vmap[v_src].R;
                emap[e].X = vmap[v_src].X;
            }
        }
    }
}

auto reset_movements(SpatialGraph &g, std::map<Vertex, vdata>& vmap, std::map<Edge, edata>& emap) -> void {
    /*
     * Send everyone back to home vertex, delete all edges and edge data.
     */
    for (int v = 0; v < g.num_vertices(); v++) {

        double summed_S = vmap[v].S * vmap[v].N;
        double summed_I = vmap[v].I * vmap[v].N;
        double summed_R = vmap[v].R * vmap[v].N;
        double summed_X = vmap[v].X * vmap[v].N;
        double summed_N = vmap[v].N;

        auto it_pair = g.out_edges(v);
        for (; it_pair.first != it_pair.second; it_pair.first++) {
            auto e = *it_pair.first; // edge
            summed_S += emap[e].S * emap[e].N;
            summed_I += emap[e].I * emap[e].N;
            summed_R += emap[e].R * emap[e].N;
            summed_X += emap[e].X * emap[e].N;
            summed_N += emap[e].N;

            //emap[e].N = 0;

        }
        vmap[v].S = summed_S/summed_N;
        vmap[v].I = summed_I/summed_N;
        vmap[v].R = summed_R/summed_N;
        vmap[v].X = summed_X/summed_N;
        vmap[v].N = summed_N;
    }
    //g.remove_all_edges();
    //emap.clear();
}

auto add_infected(SpatialGraph& g, std::map<Vertex, vdata>& vmap, double x, double y) -> void{

    Vertex v = g.closest_vertex(x, y);
    vmap[v].I = (double)(vmap[v].I*vmap[v].N+1)/vmap[v].N;
    vmap[v].S = 1-vmap[v].I;
    vmap[v].R = 0;
    vmap[v].X = 0;
    std::cout << "Added infected individual located at " << g.pos(v).first << "," << g.pos(v).second << std::endl;

}

auto gen_interactions(SpatialGraph &g, std::map<Vertex, vdata>& vmap, std::map<Edge, edata>& emap, double beta, double mu, double kappa, double alpha) -> void{
    for (int v = 0; v < g.num_vertices(); v++) {
        auto it_pair = g.in_edges(v);
        double summed_I = vmap[v].I*vmap[v].N;
        double summed_N = vmap[v].N;

        for(; it_pair.first!=it_pair.second; it_pair.first++) {
            auto e = *it_pair.first; // edge
            summed_I += emap[e].I*emap[e].N;
            summed_N += emap[e].N;
        }

        if (summed_N <= 0){
            continue;
        }
        double S_delta = -beta*vmap[v].S*summed_I/summed_N;
        double I_delta = beta*vmap[v].S*summed_I/summed_N - kappa*vmap[v].I - mu*vmap[v].I - alpha*vmap[v].I;
        double R_delta = mu*vmap[v].I + mu*vmap[v].X ;
        double D_delta = alpha*vmap[v].I + alpha*vmap[v].X;
        double X_delta = kappa*vmap[v].I - mu*vmap[v].X - alpha*vmap[v].X;

        vmap[v].S += S_delta;
        vmap[v].I += I_delta;
        vmap[v].R += R_delta;
        vmap[v].D += D_delta;
        vmap[v].X += X_delta;

        it_pair = g.in_edges(v);
        for(; it_pair.first!=it_pair.second; it_pair.first++) {
            auto e = *it_pair.first; // edge

            S_delta = -beta*emap[e].S*summed_I/summed_N;
            I_delta = beta*emap[e].S*summed_I/summed_N - kappa*emap[e].I - mu*emap[e].I - alpha*emap[e].I;
            R_delta = mu*emap[e].I + mu*emap[e].X;
            D_delta = alpha*emap[e].I + alpha*emap[e].X;
            X_delta = kappa*emap[e].I - mu*emap[e].X - alpha*emap[e].X;

            emap[e].S += S_delta;
            emap[e].I += I_delta;
            emap[e].R += R_delta;
            emap[e].D += D_delta;
            emap[e].X += X_delta;
        }
    }
}

template<typename Distribution>
auto update_waiting_times(std::map<Vertex, vdata>& vmap, Distribution dist) -> void{
    /*
     * Update the waiting times of each individual according to the random distribution dist.
     */
    for (auto& p: vmap) {
        for (auto& wt: p.second.waiting_times) {
            if (wt == 0) {
                wt = dist(global_engine())-1;
            }
            else {
                wt--;
            }
        }
    }
}

auto print_status(SpatialGraph &g, std::map<Vertex, vdata>& vmap) {
    double summed_S = 0;
    double summed_I = 0;
    double summed_R = 0;
    double summed_D = 0;
    double summed_X = 0;
    double summed_N = 0;
    for (int v = 0; v < g.num_vertices(); v++) {

        summed_S += vmap[v].S * vmap[v].N;
        summed_I += vmap[v].I * vmap[v].N;
        summed_R += vmap[v].R * vmap[v].N;
        summed_D += vmap[v].D * vmap[v].N;
        summed_X += vmap[v].X * vmap[v].N;
        summed_N += vmap[v].N;
    }
    std::cout << "S = " << summed_S << ", I = " << summed_I << ", R = " << summed_R << ", X = " << summed_X << ", D = " << summed_D << ", N = " << summed_N << "\n";
}

void to_csv(SpatialGraph& g, std::map<Vertex, vdata>& vmap, std::map<Edge, edata>, const std::string& id, std::string dir)
{
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

    for (int v = 0; v<g.num_vertices(); v++)
    {
        myfile << g.x_pos(v) << ",";
        myfile << g.y_pos(v) << ",";
        myfile << vmap[v].S*vmap[v].N <<",";
        myfile << vmap[v].I*vmap[v].N <<",";
        myfile << vmap[v].R*vmap[v].N <<",";
        myfile << vmap[v].D*vmap[v].N <<",";
        myfile << vmap[v].X*vmap[v].N <<",";
        myfile << vmap[v].N <<"\n";
    }
    myfile.close();
}

#endif //EPIGRAPH_SPATIAL_META_SIR_H
