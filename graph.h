//
// Created by roryh on 21/06/2020.
//

#ifndef EPIGRAPH_GRAPH_H
#define EPIGRAPH_GRAPH_H

#include "edge.h"

#include <iostream>
#include <fstream>
#include <random>
#include <vector>
#include <algorithm>
#include <chrono>
#include <tuple>

using Vertex = int;

class Graph {
    /*
     * A directed graph. The adjacency list is built on vectors and is optimised for quick traversal
     * and thus is slow for inserting new edges and vertices.
     *
     * Vertices are represented as unique integers.
     */
public:
    using adj_list_type = std::vector<std::vector<Edge>>;
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
    Graph() : m_num_vertices(0), m_num_edges(0) {};

    // structure modification

    auto add_vertex() -> Vertex {
        /*
         * Add a new vertex to the graph.
         */
        m_out_adj_list.emplace_back();
        m_in_adj_list.emplace_back();
        ++m_num_vertices;
        Vertex v = m_out_adj_list.size() - 1;

        return v;
    }

    auto reserve_vertices(int n) -> void {
        /*
         * Reserve memory for n vertices in the graph, this helps prevent fragmentation but is not necessary.
         */
        m_out_adj_list.reserve(n);
        m_in_adj_list.reserve(n);
    }

    virtual auto add_edge(Vertex v_src, Vertex v_dst) -> std::pair<Edge, bool> {
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
            if (m_out_adj_list.at(v_src).empty()) {
                m_out_adj_list[v_src].push_back(e);
            } else {
                m_out_adj_list[v_src].insert(out_e_it, e);
            }
            if (m_in_adj_list.at(v_dst).empty()) {
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

        for (auto &k: m_out_adj_list) {
            k.clear();
        }
        for (auto &k: m_in_adj_list) {
            k.clear();
        }


        //m_out_adj_list.resize(m_num_vertices);
        //m_in_adj_list.resize(m_num_vertices);

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
        return std::make_pair(m_in_adj_list.at(v).begin(), m_in_adj_list.at(v).end());
    }

    auto out_edge(Vertex v_src, Vertex v_dst) -> std::pair<OutEdgeIter, bool> {
        /*
         * Return an iterator that points to edge e in the out adjacency list of the source vertex of edge e. If the
         * edge does not exist then the boolean flag is false, otherwise true.
         */
        Edge e(v_src, v_dst);
        OutEdgeIter e_begin = m_out_adj_list.at(v_src).begin();
        OutEdgeIter e_end = m_out_adj_list.at(v_src).end();

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
        InEdgeIter e_begin = m_in_adj_list.at(v_dst).begin();
        InEdgeIter e_end = m_in_adj_list.at(v_dst).end();;

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

    auto num_vertices() const -> size_t {
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

    auto print_edge(const Edge &e) -> void {
        /*
         * Print the edge e -> (src, dst)
         */
        std::cout << "(" << e.src << "," << e.dst << ") ";
    }

public:

    adj_list_type m_out_adj_list; // out adjacency list
    adj_list_type m_in_adj_list; // in adjacency list
    int m_num_vertices;
    int m_num_edges;
};

template<typename VData, typename EData>
class Network {
    /*
     * A directed graph. The adjacency list is built on vectors and is optimised for quick traversal
     * and thus is slow for inserting new edges and vertices.
     *
     * Vertices are represented as unique integers.
     */
public:
    using adj_list_type = std::vector<std::vector<Edge>>;
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
    Network() : m_num_vertices(0), m_num_edges(0) {};

    // structure modification

    auto add_vertex() -> Vertex {
        /*
         * Add a new vertex to the graph.
         */
        m_out_adj_list.emplace_back();
        m_in_adj_list.emplace_back();
        ++m_num_vertices;
        Vertex v = m_out_adj_list.size() - 1;

        vprop.emplace_back(); // default construct

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
            if (m_out_adj_list.at(v_src).empty()) {
                m_out_adj_list[v_src].push_back(e);
            } else {
                m_out_adj_list[v_src].insert(out_e_it, e);
            }
            if (m_in_adj_list.at(v_dst).empty()) {
                m_in_adj_list[v_dst].push_back(e);
            } else {
                m_in_adj_list[v_dst].insert(in_e_it, e);
            }

            // add edge property
            eprop[e]; // default construct
            return std::make_pair(e, true);
        }
        else {
            return std::make_pair(e, false);
        }
    }

    auto remove_all_edges() -> void {
        /*
        Removes all edges from the graph.
         */

        for (auto &k: m_out_adj_list) {
            k.clear();
        }
        for (auto &k: m_in_adj_list) {
            k.clear();
        }


        //m_out_adj_list.resize(m_num_vertices);
        //m_in_adj_list.resize(m_num_vertices);

        m_num_edges = 0;

        eprop.clear();
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
        return std::make_pair(m_in_adj_list.at(v).begin(), m_in_adj_list.at(v).end());
    }

    auto out_edge(Vertex v_src, Vertex v_dst) -> std::pair<OutEdgeIter, bool> {
        /*
         * Return an iterator that points to edge e in the out adjacency list of the source vertex of edge e. If the
         * edge does not exist then the boolean flag is false, otherwise true.
         */
        Edge e(v_src, v_dst);
        OutEdgeIter e_begin = m_out_adj_list.at(v_src).begin();
        OutEdgeIter e_end = m_out_adj_list.at(v_src).end();

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
        InEdgeIter e_begin = m_in_adj_list.at(v_dst).begin();
        InEdgeIter e_end = m_in_adj_list.at(v_dst).end();;

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

    auto num_vertices() const -> size_t {
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

    auto print_edge(const Edge &e) -> void {
        /*
         * Print the edge e -> (src, dst)
         */
        std::cout << "(" << e.src << "," << e.dst << ") ";
    }

public:

    adj_list_type m_out_adj_list; // out adjacency list
    adj_list_type m_in_adj_list; // in adjacency list
    int m_num_vertices;
    int m_num_edges;

    std::vector<VData> vprop;
    std::map<Edge, EData> eprop;
};

#endif //EPIGRAPH_GRAPH_H
