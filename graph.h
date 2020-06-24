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

    auto add_vertex() -> Vertex;

    auto reserve_vertices(int n) -> void;

    auto add_edge(Vertex v_src, Vertex v_dst) -> std::pair<Edge, bool>;

    auto remove_all_edges() -> void;

    // structure queries
    auto out_edges(Vertex v) -> std::pair<OutEdgeIter, OutEdgeIter>;

    auto in_edges(Vertex v) -> std::pair<InEdgeIter, InEdgeIter>;

    auto out_edge(Vertex v_src, Vertex v_dst) -> std::pair<OutEdgeIter, bool>;

    auto in_edge(Vertex v_src, Vertex v_dst) -> std::pair<InEdgeIter, bool>;

    auto vertex(Vertex v) -> bool;

    auto edge(Vertex v_src, Vertex v_dst) -> std::pair<Edge, bool>;

    auto num_vertices() const -> size_t;

    auto num_edges() -> size_t;

    auto print_edge(const Edge &e) -> void;


public:

    adj_list_type m_out_adj_list; // out adjacency list
    adj_list_type m_in_adj_list; // in adjacency list
    int m_num_vertices;
    int m_num_edges;
};


#endif //EPIGRAPH_GRAPH_H
