// list.h
// Declaration of the class List.

#ifndef LIST_H
#define LIST_H

#include <iostream>
#include <string>
#include "allele.h"

class List
{
  friend ostream & operator<<( ostream & out, const List & list );

  public:
    List();
    ~List();
    List & insert( Allele * ap );                    // add to end of list
    List & operator=( const List & list2 );
    Allele * search( string name );
    Allele * search_by_conversion( int x );

    Allele * get_first() const { return first; }
    Allele * get_last() const { return last; }

  private:
    Allele * first;                       // pointer to first allele in list
    Allele * last;                        // pointer to last allele in list
};

#endif
