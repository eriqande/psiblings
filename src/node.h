// node.h
// Declaration of the class Node.

#ifndef NODE_H
#define NODE_H

#include <iostream>
#include <string>
using namespace std;

class Edge;                    // forward declaration

class Node
{
  public:
    Node();
    Node( const Node & x );
    ~Node();
    Node & operator=( const Node & x );
                    
    int get_id() const { return id; }
    string get_trueID() const { return trueID; }
    Edge * get_first() const { return first_edge; }
    Edge * get_scan() const { return scan_ptr; }
    Node * get_parent() const { return parent; }
    int get_mincap() const { return mincap; }
    int get_subtreeSize() const { return subtreeSize; }
    bool get_memberT() const { return memberT; }
    bool get_checked() const { return checked; }
    int get_dist() const { return dist; }
    int get_excess() const { return excess; }
    Node * get_bfs_link() const { return bfs_link; }
    Node * get_stack_link() const { return stack_link; }
    bool get_alive() const { return alive; }
    bool get_unmarked() const { return unmarked; }

    void set_id( int x ) { id = x; }
    void set_trueID( string s ) { trueID = s; }
    void set_first( Edge * ep ) { first_edge = ep; }
    void set_scan( Edge * ep ) { scan_ptr = ep; }
    void set_parent( Node * np ) { parent = np; }
    void set_mincap( int x ) { mincap = x; }
    void set_subtreeSize( int x ) { subtreeSize = x; }
    void set_memberT( bool x ) { memberT = x; }
    void set_checked( bool x ) { checked = x; }
    void set_dist( int x ) { dist = x; }
    void set_excess( int x ) { excess = x; }
    void set_bfs_link( Node * np ) { bfs_link = np; }
    void set_stack_link( Node * np ) { stack_link = np; }
    void set_alive( bool x ) { alive = x; }
    void set_unmarked( bool x ) { unmarked = x; }
       
  private:
    int id;                    // id number for use in ghc-tree
    string trueID;             // id given in user data file
    Edge * first_edge;         // first in list of incident edges
    Edge * scan_ptr;           // next to be scanned when node visited again
    Node * parent;             // pointer for Gomory-Hu cut tree
    int mincap;                // capacity of min cut between node and its parent
    int subtreeSize;           // size of subtree including this node
    bool memberT;              // true if node is on T side of min cut
    bool checked;              // true if mincut score already examined

    // below for use by maxflow()
    int dist;
    int excess;
    Node * bfs_link;           // for bfs queue
    Node * stack_link;         // for stack of active nodes
    bool alive;
    bool unmarked;             // while bfs in progress
};

#endif
