#include <R.h>
// edge.cc
// Function definitions for class Edge.


#include "edge.h"


/*-----------------------------  constructor  -----------------------------*/

Edge::Edge() {
  adjacent = 0;
  next = 0;
  back = 0;
  cap = 0;
  rcap = 0;
}



/*--------------------------  copy constructor  ---------------------------*/

Edge::Edge( const Edge & x ) {
  adjacent = x.adjacent;
  next = x.next;
  back = x.back;
  cap = x.cap;
  rcap = x.rcap;
}



/*------------------------------  destructor  -----------------------------*/

Edge::~Edge() {
  adjacent = 0;
  next = 0;
  back = 0;
  cap = 0;
  rcap = 0;
}



/*----------------------------  assignment  -------------------------------*/

Edge & Edge::operator=( const Edge & x ) {
  if( this == &x ) return *this;

  adjacent = x.adjacent;
  next = x.next;
  back = x.back;
  cap = x.cap;
  rcap = x.rcap;

  return *this;
}
