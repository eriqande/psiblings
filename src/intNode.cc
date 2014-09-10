#include <R.h>
// intNode.cc
// Function definitions for class IntNode.


#include "intNode.h"


/*-----------------------------  constructor  -----------------------------*/

IntNode::IntNode( int id ) {
  id_num = id;
  notInList = false;
  next = 0;
}



/*--------------------------  copy constructor  ---------------------------*/

IntNode::IntNode( const IntNode & x ) {
  id_num = x.id_num;
  notInList = x.notInList;
  next = x.next;
}



/*------------------------------  destructor  -----------------------------*/

IntNode::~IntNode() {
  id_num = 0;
  notInList = false;
  next = 0;
}



/*----------------------------  assignment  -------------------------------*/

IntNode & IntNode::operator=( const IntNode & x ) {
  if( this == &x ) return *this;

  id_num = x.id_num;
  notInList = x.notInList;
  next = x.next;

  return *this;
}



/*-----------------------------  operator==  ------------------------------*/

bool IntNode::operator==( const IntNode & x ) const {
  if( id_num == x.id_num )
    return true;
  else return false;
}



/*--------------------------------  print()  ------------------------------*/

void IntNode::print( ostream & out ) const {
  out << id_num << " ";
}



/*--------------------------------  copy()  -------------------------------*/

IntNode * IntNode::copy() const {
  IntNode * ip;
  ip = new IntNode( id_num );

  ip->set_next( next );

  return ip;
}
