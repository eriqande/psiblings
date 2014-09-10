#include <R.h>
// list.cc
// Function definitions for class List.


#include "list.h"


/*----------------------------  constructor  -----------------------------*/

List::List() {
  first = 0;
  last = 0;
}



/*-----------------------------  destructor  -----------------------------*/

List::~List() {
  if( first != 0 )  {                       // if list is not empty
    Allele * currentPtr;
    Allele * tempPtr;
    currentPtr = first;
    tempPtr = currentPtr;
    last = 0;

    while( tempPtr != 0 ) {                 // delete remaining nodes
      currentPtr = tempPtr->get_next();
      delete tempPtr;
      tempPtr = currentPtr;
    }
  }
}



/*-------------------------------  insert()  -------------------------------
|   Takes an allele pointer as a parameter.                                 |
|   Allele pointer ap must point to an allele created with 'new'.           |
|   Returns a list with the new allele at end of list.                      |
--------------------------------------------------------------------------*/

List & List::insert(Allele * ap) {
  if( first == 0 ) {                        // if list is empty
    ap->set_next(0);
    first = ap;
    last = ap;
    return *this;
  }

  else {                                    // list is not empty
    ap->set_next( 0 );
    last->set_next( ap );
    last = ap;
  }

  return *this;
}



/*-----------------------------  assignment  -------------------------------
|   Assigns the object list2 to the calling object.                        |
|   A reference to the calling object is returned to enable cascading.     |
--------------------------------------------------------------------------*/

List & List::operator=( const List & list2 ) {
  Allele * currentPtr;
  Allele * tempPtr;
  Allele * mp;
  Allele * mp2;

  if( &list2 != this ) {                    // check for self-assignment
    if( first ) {                           // if list is not empty
      currentPtr = first;

      while( currentPtr ) {                 // prevent a memory leak
        tempPtr = currentPtr;               // delete alleles in calling list
        currentPtr = currentPtr->get_next();
        delete tempPtr;
      }
    }

    //set calling list's first pointer to a copy of first allele in list2
    currentPtr = list2.first;
    mp = currentPtr->copy();                // copy first allele in list2
    first = mp;                             // set first to the new copy
    currentPtr = currentPtr->get_next();    // advance currentPtr

    //copy rest of list2 and link alleles together
    while( currentPtr ) {
      mp2 = currentPtr->copy();             // copy next allele in list2
      mp->set_next( mp2 );                  // set prev Next to new allele
      last = mp2;                           // set last pointer
      mp = mp2;                             // advance mp
      currentPtr = currentPtr->get_next();  // advance currentPtr
    }
  }

  return *this;
}



/*------------------------------  search()  --------------------------------
|   Searches a list for an allele specified by name.                       |
|   A pointer to the allele or a NULL pointer is returned.                 |
--------------------------------------------------------------------------*/

Allele * List::search( string name ) {
  Allele * currentPtr;
  string s1;

  // traverse the list looking for a match to the specified name
  for(currentPtr = first; currentPtr; currentPtr = currentPtr->get_next()) {
    s1 = currentPtr->get_name();
    if( s1.compare( name ) == 0 )        // if a match is found
      return currentPtr;                 // return a pointer to that allele
  }

  // no match found
  return currentPtr;                     // currentPtr is null(end of list)
}



/*------------------  search_by_conversion()  ------------------------------
|   Searches a list for an allele specified by its conversion number.      |
|   A pointer to the allele or a NULL pointer is returned.                 |
--------------------------------------------------------------------------*/

Allele * List::search_by_conversion( int x ) {
  Allele * currentPtr;
  int y;

  // traverse the list looking for a match to the specified name
  for(currentPtr = first; currentPtr; currentPtr = currentPtr->get_next()) {
    y = currentPtr->get_conversion();
    if( y == x )                         // if a match is found
      return currentPtr;                 // return a pointer to that allele
  }

  // no match found
  return currentPtr;                     // currentPtr is null(end of list)
}



/*-----------------------------  operator<<  -------------------------------
|   Prints the alleles in the list.                                        |
|   Used for debugging.                                                    |
--------------------------------------------------------------------------*/

ostream & operator<<( ostream & out, const List & list ) {
  Allele * currentPtr = list.first;

  if( list.first == 0 ) {
    out << "The list is empty." << endl;
    return out;
  }

  while( currentPtr != 0 ) {
    currentPtr->print( out );
    currentPtr = currentPtr->get_next();
  }
  
  return out;
}
