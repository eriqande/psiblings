// graph.h
// Declaration of the class Graph.

#ifndef GRAPH_H
#define GRAPH_H

#include <iostream>
#include "node.h"
#include "edge.h"

class Graph
{
  public:
    Graph( int n, int m );
    Graph( const Graph & x );
    ~Graph();
    Graph & operator=( const Graph & x );
                     
    int get_n_nodes() const { return n_nodes; }
    int get_n_edges() const { return n_edges; }
    int get_n_edges0() const { return n_edges0; }
    Node * get_nodes() const { return nodes; }
    Edge * get_edges() const { return edges; }

    void set_n_nodes( int x ) { n_nodes = x; }
    void set_n_edges( int x ) { n_edges = x; }
    void set_n_edges0( int x ) { n_edges0 = x; }
    void set_nodes( Node * np ) { nodes = np; }
    void set_edges( Edge * ep ) { edges = ep; }

    Node * nodes;                    // pointer to array of Nodes
    Edge * edges;                    // pointer to list of Edges

  private:
    int n_nodes;                     // number of nodes
    int n_edges;                     // number of non-zero capacity edges
    int n_edges0;                    // total number of edges
};

#endif
