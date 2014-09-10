#include <R.h>
// graph.cc
// Function definitions for class Graph.


#include "graph.h"


/*-----------------------------  constructor  -----------------------------*/

Graph::Graph( int n, int m ) {
  n_nodes = n;
  n_edges = m;
  n_edges0 = m;
  nodes = new Node[n];
  edges = new Edge[2 * m];
}



/*--------------------------  copy constructor  ---------------------------*/

Graph::Graph( const Graph & x ) {
  n_nodes = x.n_nodes;
  n_edges = x.n_edges;
  n_edges0 = x.n_edges0;
  nodes = x.nodes;
  edges = x.edges;
}



/*------------------------------  destructor  -----------------------------*/

Graph::~Graph() {
  n_nodes = 0;
  n_edges = 0;
  n_edges0 = 0;
  nodes = 0;
  edges = 0;
}



/*----------------------------  assignment  -------------------------------*/

Graph & Graph::operator=( const Graph & x ) {
  if( this == &x ) return *this;                 // check for self assignment

  n_nodes = x.n_nodes;
  n_edges = x.n_edges;
  n_edges0 = x.n_edges0;
  nodes = x.nodes;
  edges = x.edges;

  return *this;
}
