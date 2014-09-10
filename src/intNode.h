// intNode.h
// Declaration of the class IntNode.

#ifndef INTNODE_H
#define INTNODE_H

#include <iostream>
#include <string>

using namespace std;

class IntNode
{
  public:
    IntNode( int id );
    IntNode( const IntNode & x );
    ~IntNode();

    IntNode & operator=( const IntNode & x );
    bool operator==( const IntNode & x ) const;

    void print( ostream & out ) const;
    IntNode * copy() const;

    int get_id_num() const { return id_num; }
    int get_notInList() const { return notInList; }
    IntNode * get_next() const { return next; }

    void set_notInList( bool x ) { notInList = x; }
    void set_next( IntNode * ip ) { next = ip; }

  private:
    int id_num;                        // integer alias for individual name
    bool notInList;                    // used to temporarily remove from list
    IntNode * next;
};

#endif
