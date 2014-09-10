// intList.h
// Declaration of the class IntList.

#ifndef INT_LIST_H
#define INT_LIST_H

#include <iostream>
#include <string>
#include "intNode.h"


using namespace std; 

class intList
{
  friend ostream & operator<<( ostream & out, const intList & list );

  public:
    intList();
    ~intList();
    intList & insert( IntNode * ip );            // add to end of list
    intList & operator=( const intList & list2 );

    bool isEmpty() const;
    void mark_STcut( int set, int cc[] );
    intList & delete_first();
    void delete_node( IntNode * ip );
    void delete_notInList();
    void restore_notInList();
    int countAdj( int x );                       // count neighbors with id > x
    bool search( int  x );
    IntNode * get_first() const { return first; }
    IntNode * get_last() const { return last; }

  private:
    IntNode * first;                          // pointer to first allele in list
    IntNode * last;                           // pointer to last allele in list
};

#endif
