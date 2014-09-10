// edge.h
// Declaration of the class Edge.

#ifndef EDGE_H
#define EDGE_H

#include <iostream>

class Node;                          // forward declaration

class Edge
{
  public:
    Edge();
    Edge( const Edge & x );
    ~Edge();
    Edge & operator=( const Edge & x );
                     
    Node * get_adjacent() const { return adjacent; }
    Edge * get_next() const { return next; }
    Edge * get_back() const { return back; }
    int get_cap() const { return cap; }
    int get_rcap() const { return rcap; }

    void set_adjacent( Node * np ) { adjacent = np; }
    void set_next( Edge * ep ) { next = ep; }
    void set_back( Edge * ep ) { back = ep; }
    void set_cap( int x ) { cap = x; }
    void set_rcap( int x ) { rcap = x; }

  private:
    Node * adjacent;                 // pointer to adjacent node
    Edge * next;                     // next edge in incidence list for node
    Edge * back;                     // pointer to reverse edge
    int cap;                         // capacity
    int rcap;                        // residual capacity
};

#endif
