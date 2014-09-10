#include <R.h>
// node.cc
// Function definitions for class Node.


#include "node.h"


/*-----------------------------  constructor  -----------------------------*/

Node::Node() {
  id = 0;
  trueID = "";
  first_edge = scan_ptr = 0;
  parent = 0;
  mincap = 0;
  subtreeSize = 1;
  memberT = checked = false;
   
  dist = excess = 0;
  bfs_link = stack_link = 0;
  alive = unmarked = 0;
}



/*--------------------------  copy constructor  ---------------------------*/

Node::Node( const Node & x ) {
  id = x.id;
  trueID = x.trueID;
  first_edge = x.first_edge;
  scan_ptr = x.scan_ptr;
  parent = x.parent;
  mincap = x.mincap;
  subtreeSize = x.subtreeSize;
  memberT = x.memberT;
  checked = x.checked;

  dist = x.dist;
  excess = x.excess;
  bfs_link = x.bfs_link;
  stack_link = x.stack_link;
  alive = x.alive;
  unmarked = x.unmarked;
}



/*------------------------------  destructor  -----------------------------*/

Node::~Node() {
  id = 0;
  trueID = "";
  first_edge = scan_ptr = 0;
  parent = 0;
  mincap = subtreeSize = 0;
  memberT = checked = false;

  dist = excess = 0;
  bfs_link = stack_link = 0;
  alive = unmarked = 0;
}



/*----------------------------  assignment  -------------------------------*/

Node & Node::operator=( const Node & x ) {
  if( this == &x ) return *this;

  id = x.id;
  trueID = x.trueID;
  first_edge = x.first_edge;
  scan_ptr = x.scan_ptr;
  parent = x.parent;
  mincap = x.mincap;
  subtreeSize = x.subtreeSize;
  memberT = x.memberT;
  checked = x.checked;

  dist = x.dist;
  excess = x.excess;
  bfs_link = x.bfs_link;
  stack_link = x.stack_link;
  alive = x.alive;
  unmarked = x.unmarked;

  return *this;
}
