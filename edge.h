//
// Created by roryh on 05/04/2020.
//

#ifndef EPIGRAPH_EDGE_H
#define EPIGRAPH_EDGE_H

using Vertex = int;
class Edge {
    /*
     * Models a directed edge src -> dst.
     */
public:
    Vertex src;
    Vertex dst;

    Edge() {}
    Edge(Vertex v_src, Vertex v_dst): src(v_src),  dst(v_dst) {}
    Edge(const Edge& e) {
        src = e.src;
        dst = e.dst;
    }
    void operator=(Edge e){
        src = e.src;
        dst = e.dst;
    }
    bool operator<(const Edge& e) const{
        return src < e.src || (src == e.src && dst < e.dst);
    }
    bool operator==(const Edge& e) {
        return src == e.src && dst == e.dst;
    }
    bool operator!=(const Edge &e) {
        return !(*this == e);
    }

    //int id() {return m_id;}
};

#endif //EPIGRAPH_EDGE_H
