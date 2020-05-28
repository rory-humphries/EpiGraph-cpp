//
// Created by roryh on 05/04/2020.
//

#ifndef EPIGRAPH_SPATIAL_GRAPH_H
#define EPIGRAPH_SPATIAL_GRAPH_H

#include "distributions.h"
#include "edge.h"


#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <iostream>
#include <fstream>
#include <random>
#include <vector>
#include <algorithm>
#include <chrono>
const double pi = 3.14;


namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;
using Vertex = int;

class SpatialGraph {
    /*
     * A directed graph where each vertex has a spatial (x,y) index stored in a R*-tree data structure. The adjacency
     * list is built on vectors and is optimised for quick traversal and thus is slow for inserting new edges and
     * vertices.
     *
     * Vertices are represented as unique integers.
     */
public:
    // iterators
    using OutEdgeIter = std::vector<Edge>::iterator;
    using InEdgeIter = std::vector<Edge>::iterator;

    class EdgeIter {
        /*
         * An edge iterator class which iterates over the entire edge set.
         */
    public:
        using difference_type = ptrdiff_t;
        using value_type = Edge;
        using pointer = Edge *;
        using reference = Edge &;
        using iterator_category = std::forward_iterator_tag;

        EdgeIter(std::vector<Edge>::iterator e_curr,
                 std::vector<Edge>::iterator e_end,
                 std::vector<std::vector<Edge>>::iterator v_curr,
                 std::vector<std::vector<Edge>>::iterator v_end) : m_e_curr(e_curr), m_e_end(e_end),
                                                                   m_edges(v_curr),
                                                                   m_edges_end(v_end) {
            if (m_edges->empty()) {
                while (m_edges != m_edges_end && m_edges->empty()) {
                    m_edges++;
                }
                m_e_curr = m_edges->begin();
                m_e_end = m_edges->end();
            }
        }

        EdgeIter &operator++() {
            m_e_curr++;
            if (m_e_curr != m_e_end) {
                return *this;
            } else {
                m_edges++;
                while (m_edges != m_edges_end && m_edges->empty()) {
                    m_edges++;
                }
                m_e_curr = m_edges->begin();
                m_e_end = m_edges->end();
                return *this;
            }
        }

        EdgeIter operator++(int) {
            EdgeIter retval = *this;
            ++(*this);
            return retval;
        }

        bool operator==(EdgeIter other) const {
            return m_e_curr == other.m_e_curr && m_edges == other.m_edges;
        }

        bool operator!=(EdgeIter other) const {
            return !(*this == other);
        }

        reference operator*() {
            return *(m_e_curr);
        }
        // iterator traits
    private:
        std::vector<Edge>::iterator m_e_end;
        std::vector<Edge>::iterator m_e_curr;
        std::vector<std::vector<Edge>>::iterator m_edges;
        std::vector<std::vector<Edge>>::iterator m_edges_end;
    };

    // constructors

    // Default constructor
    SpatialGraph() : m_num_vertices(0), m_num_edges(0), rtree() {};

    // structure modification

    auto add_vertex(double x, double y) -> Vertex {
        /*
         * Add a new vertex to the graph at position (x,y).
         */
        m_out_adj_list.emplace_back();
        m_in_adj_list.emplace_back();
        ++m_num_vertices;
        Vertex v = m_out_adj_list.size() - 1;

        point p(x, y);
        m_positions.push_back(p);
        rtree.insert(std::make_pair(p, v)); // use vertex index
        return v;
    }

    auto reserve_vertices(int n) -> void {
        /*
         * Reserve memory for n vertices in the graph, this helps prevent fragmentation but is not necessary.
         */
        m_out_adj_list.reserve(n);
        m_in_adj_list.reserve(n);
    }

    auto add_edge(Vertex v_src, Vertex v_dst) -> std::pair<Edge, bool> {
        /*
         * Add a directed edge connecting v_src to v_dst. Returns a pair <Edge, bool>, if the boolean variable is true
         * the edge was added, if false, then it already exists.
         */
        Edge e(v_src, v_dst);
        bool edge_added = false;

        OutEdgeIter out_e_it;
        bool out_edge_exists;
        std::tie(out_e_it, out_edge_exists) = out_edge(v_src, v_dst);

        InEdgeIter in_e_it;
        bool in_edge_exists;
        std::tie(in_e_it, in_edge_exists) = in_edge(v_src, v_dst);

        if (!out_edge_exists && !in_edge_exists) {
            if (m_out_adj_list[v_src].empty()) {
                m_out_adj_list[v_src].push_back(e);
            } else {
                m_out_adj_list[v_src].insert(out_e_it, e);
            }
            if (m_in_adj_list[v_dst].empty()) {
                m_in_adj_list[v_dst].push_back(e);
            } else {
                m_in_adj_list[v_dst].insert(in_e_it, e);
            }
            return std::make_pair(e, true);
        }
        return std::make_pair(e, false);
    }

    auto remove_all_edges() -> void {
        /*
        Removes all edges from the graph.
         */
        m_out_adj_list.clear();
        m_in_adj_list.clear();

        m_out_adj_list.resize(m_num_vertices);
        m_in_adj_list.resize(m_num_vertices);

        m_num_edges = 0;
    }

    // structure queries
    auto out_edges(Vertex v) -> std::pair<OutEdgeIter, OutEdgeIter> {
        /*
         * Return an iterator that traverses all out edges from vertex v.
         */
        return std::make_pair(m_out_adj_list[v].begin(), m_out_adj_list[v].end());
    }

    auto in_edges(Vertex v) -> std::pair<InEdgeIter, InEdgeIter> {
        /*
         * Return an iterator that traverses all in edges to vertex v.
         */
        return std::make_pair(m_in_adj_list[v].begin(), m_in_adj_list[v].end());
    }

    auto out_edge(Vertex v_src, Vertex v_dst) -> std::pair<OutEdgeIter, bool> {
        /*
         * Return an iterator that points to edge e in the out adjacency list of the source vertex of edge e. If the
         * edge does not exist then the boolean flag is false, otherwise true.
         */
        Edge e(v_src, v_dst);
        OutEdgeIter e_begin = m_out_adj_list[v_src].begin();
        OutEdgeIter e_end = m_out_adj_list[v_src].end();

        // find first edge
        auto e_pos = std::lower_bound(e_begin, e_end, e);
        if (std::distance(e_pos, e_end) == 0) {
            return std::make_pair(e_pos, false);
        } else {
            return std::make_pair(e_pos, *e_pos == e);
        }
    }

    auto in_edge(Vertex v_src, Vertex v_dst) -> std::pair<InEdgeIter, bool> {
        /*
        * Return an iterator that points to edge e in the out adjacency list of the source vertex of edge e. If the
        * edge does not exist then the boolean flag is false, otherwise true.
        */
        Edge e(v_src, v_dst);
        InEdgeIter e_begin = m_in_adj_list[v_dst].begin();
        InEdgeIter e_end = m_in_adj_list[v_dst].end();;

        // find first edge
        auto e_pos = std::lower_bound(e_begin, e_end, e);
        if (std::distance(e_pos, e_end) == 0) {
            return std::make_pair(e_pos, false);
        } else {
            return std::make_pair(e_pos, *e_pos == e);
        }
    }

    auto vertex(Vertex v) -> bool {
        /*
         * Return true if vertex v is in the graph.
         */
        if (v < num_vertices() && v > 0)
            return true;
        else
            return false;
    }

    auto edge(Vertex v_src, Vertex v_dst) -> std::pair<Edge, bool> {
        /*
         * Return the edge that connects v_src to v_dst, if the boolean flag is true the edge exists.
         */
        auto e_begin = m_out_adj_list[v_src].begin();
        auto e_end = m_out_adj_list[v_src].end();

        // find first edge
        Edge e(v_src, v_dst);
        auto e_pos = std::lower_bound(e_begin, e_end, e);
        // add a check
        return std::make_pair(e, *e_pos == e);
    }

    auto num_vertices() -> size_t {
        /*
         * Return the number of vertices in the graph.
         */
        return m_num_vertices;
    }

    auto num_edges() -> size_t {
        /*
         * Return the number of edges in the graph.
         */
        return m_num_edges;
    }

    // spacial queries
    auto distance(Vertex v1, Vertex v2) -> double {
        /*
         * Return the distance between two vertices v1 and v2.
         */
        return bg::distance(m_positions[v1], m_positions[v2]);
    }
    auto long_lat_distance(Vertex v1, Vertex v2) -> double {

        double lat1 = m_positions[v1].y();
        double lon1 = m_positions[v1].x();
        double lat2 = m_positions[v2].y();
        double lon2 = m_positions[v2].x();

        double R = 6371e3; // metres
        double phi1 = lat1 * pi / 180.0; // φ, λ in radians
        double phi2 = lat2 * pi / 180.0;
        double delphi = (lat2 - lat1) * pi / 180.0;
        double dellam = (lon2 - lon1) * pi / 180.0;

        double a = sin(delphi / 2) * sin(delphi / 2) +
                  cos(phi1) * cos(phi2) *
                  sin(dellam / 2) * sin(dellam / 2);
        double c = 2 * atan2(sqrt(a), sqrt(1 - a));

        double d = R * c; // in metres
        return d;
    }
    template<typename OutputIter>
    auto vertices_at_dist(Vertex v, double radius, double err, OutputIter result) -> OutputIter {
        /*
         * Return all vertices between the distances of radius - err and radius + r.
         */

        typedef bg::model::polygon<point> ring_t;

        std::vector<std::pair<int, int>> circle;
        ring_t outer_ring;
        ring_t inner_ring;


        // construct outer ring
        bresenhams_circle(radius + err, x_pos(v), y_pos(v), std::back_inserter(circle));

        for (const auto &p: circle) {
            bg::append(outer_ring, point(p.first, p.second));
        }
        bg::correct(outer_ring);

        // construct inner ring
        // clear the circle;
        circle.clear();

        bresenhams_circle(radius - err, x_pos(v), y_pos(v), std::back_inserter(circle));
        outer_ring.inners().resize(1);
        for (const auto &p: circle) {
            //bg::append(inner_ring, point(p.first, p.second));
            bg::append(outer_ring.inners()[0], point(p.first, p.second));
        }
        bg::correct(outer_ring);

        std::vector<rtvalue> results;
        rtree.query(bgi::covered_by(outer_ring), std::back_inserter(results));
        for (auto &r: results) {
            *result = r.second;
            result++;
        }

        return result;
    }

    template<typename OutputIter>
    auto vertices_within_dist(Vertex v, double radius, OutputIter result) -> OutputIter {
        /*
         * Return all vertices within a distance of radius +- err.
         */

        typedef bg::model::polygon<point> poly;

        std::vector<rtvalue> results;
        rtree.query(bgi::satisfies(
                [&](rtvalue const &v_dst) { return bg::distance(m_positions[v], v_dst.first) < radius; }),
                    std::back_inserter(results));

        for (auto &r: results) {
            *result = r.second;
            result++;
        }

        return result;
    }

    auto closest_vertex(double x, double y) -> Vertex {
        /*
         * Find the closest vertex to the point (x,y).
         */

        rtvalue v; // pair<point, vertex>
        rtvalue *v_ptr = &v; // rtree query must take an output iterator
        rtree.query(bgi::nearest(point(x, y), 1), v_ptr);

        return v.second;
    }

    auto print_edge(const Edge &e) -> void {
        /*
         * Print the edge e -> (src, dst)
         */
        std::cout << "(" << e.src << "," << e.dst << ") ";
    }

    auto x_pos(Vertex v) -> double {
        /*
         * Return the x position of vertex v.
         */
        return m_positions[v].x();
    }

    auto y_pos(Vertex v) -> double {
        /*
         * Return the y position of vertex v.
         */
        return m_positions[v].y();
    }

    auto pos(Vertex v) -> std::pair<double, double> {
        /*
         * Return the (x,y) position of vertex v.
         */
        auto p = m_positions[v];
        return std::make_pair(p.x(), p.y());
    }

public:

    using point = bg::model::d2::point_xy<double>;// mdoels a cartesean (x,y) point
    using rtvalue = std::pair<point, Vertex>;

    std::vector<std::vector<Edge>> m_out_adj_list; // out adjacency list
    std::vector<std::vector<Edge>> m_in_adj_list; // in adjacency list
    std::vector<point> m_positions; // the positions of each vertex
    int m_num_vertices;
    int m_num_edges;
    bgi::rtree<rtvalue, bgi::quadratic<16> > rtree; // spatial index for vertices (reverse lookup to m_positions)

};

template<typename Distribution>
void build_edges(SpatialGraph &g, Distribution dist) {
    /*
     * Populate the graph with random edges whose distance are determined by the random distribution dist.
     * Dist must provide operator()(gen) where gen is a random number generator, such as meresenne_twister
     */

    std::uniform_real_distribution<> uni_dist(0, 1);

    for (int v = 0; v < g.num_vertices(); v++) {
        double d = dist(global_engine());
        double u = uni_dist(global_engine());

        Vertex vdst = g.closest_vertex(g.x_pos(v) + d * cos(2 * 3.14 * u), g.y_pos(v) + d * sin(2 * 3.14 * u));
        g.add_edge(v, vdst);
    }
}

#endif //EPIGRAPH_SPATIAL_GRAPH_H
