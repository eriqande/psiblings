#include <R.h>
// intList.cc
// Function definitions for class IntList.


#include "intList.h"

using namespace std;


/*----------------------------  constructor  -----------------------------*/

intList::intList() {
  first = 0;
  last = 0;
}



/*-----------------------------  destructor  -----------------------------*/

intList::~intList() {
  if( first != 0 )  {                         // if list is not empty
    IntNode * currentPtr;
    IntNode * tempPtr;
    currentPtr = first;
    tempPtr = currentPtr;
    last = 0;

    while( tempPtr != 0 ) {                   // delete remaining nodes
      currentPtr = tempPtr->get_next();
      delete tempPtr;
      tempPtr = currentPtr;
    }
  }
}



/*-------------------------------  insert()  -------------------------------
|   Takes an allele pointer as a parameter.                                |
|   IntNode pointer ip must point to a node created with 'new'.            |
|   Returns a list with the new node at end of list.                       |
--------------------------------------------------------------------------*/

intList & intList::insert( IntNode * ip ) {

  if( first == 0 ) {                          // if list is empty
    ip->set_next( 0 );
    first = ip;
    last = ip;
    return *this;
  }

  else {                                      // list is not empty
    ip->set_next( 0 );
    last->set_next( ip );
    last = ip;
  }

  return *this;
}



/*-----------------------------  assignment  -------------------------------
|   Assigns the object list2 to the calling object.                        |
|   A reference to the calling object is returned to enable cascading.     |
--------------------------------------------------------------------------*/

intList & intList::operator=( const intList & list2 ) {
  IntNode * currentPtr;
  IntNode * tempPtr;
  IntNode * mp;
  IntNode * mp2;

  if( &list2 != this ) {                      // check for self-assignment
    if( first ) {                             // if calling list is not empty
      currentPtr = first;

      while( currentPtr ) {                   // prevent a memory leak
        tempPtr = currentPtr;               // delete alleles in calling list
        currentPtr = currentPtr->get_next();
        delete tempPtr;
      }
    }

    //set calling list's first pointer to a copy of first allele in list2
    currentPtr = list2.first;
    if( currentPtr ) {
      mp = currentPtr->copy();                // copy first allele in list2
      first = mp;                             // set first to the new copy
      currentPtr = currentPtr->get_next();    // advance currentPtr
    }
    else {                                    // list2 is empty
      first = NULL;
      return *this;
    }

    //copy rest of list2 and link alleles together
    while( currentPtr ) {
      mp2 = currentPtr->copy();               // copy next allele in list2
      mp->set_next( mp2 );                    // set prev Next to new allele
      last = mp2;                             // set last pointer
      mp = mp2;                               // advance mp
      currentPtr = currentPtr->get_next();    // advance currentPtr
    }
  }

  return *this;
}



/*-----------------------------  isEmpty()  --------------------------------
|   Returns true if the intList is empty.                                  |
--------------------------------------------------------------------------*/

bool intList::isEmpty() const {
  if( first == NULL )
    return true;
  else
    return false;
}



/*---------------------------  mark_STcut()  -------------------------------
|   Mark nodes as 'not in list' if they are not in set 'set' (S or T).     |
--------------------------------------------------------------------------*/

void intList::mark_STcut( int set, int cc[] ) {
  int ident;
  IntNode * ip;

  for( ip = first; ip; ip = ip->get_next() ) {
    ident = ip->get_id_num();
    if( cc[ident] != set )
      ip->set_notInList( true );
  }
}



/*------------------------  delete_first()  --------------------------------
|   Removes the first element in the list.                                 |
--------------------------------------------------------------------------*/

intList & intList::delete_first() {
  IntNode * ip;

  ip = first;
  first = ip->get_next();
  if( first == NULL )
    last = NULL;
  delete ip;

  return *this;
}



/*------------------------  delete_node()  ---------------------------------
|   Deletes the node pointed to by ip.                                     |
|   Called when ip != first.                                               |
--------------------------------------------------------------------------*/

void intList::delete_node( IntNode * ip ) {
  IntNode * prev_ptr;

  for( prev_ptr = first; prev_ptr; prev_ptr = prev_ptr->get_next() ) {
    if( prev_ptr->get_next() == ip ) {
      prev_ptr->set_next( ip->get_next() );
      if( ip == last )
        last = prev_ptr;
      delete ip;
      return;
    }
  }
}



/*-------------------------  delete_notInList()  ---------------------------
|  Deletes all nodes in list that have data member 'notInList' == true.    |
--------------------------------------------------------------------------*/

void intList::delete_notInList() {
  IntNode * ip, * temp;

  ip = first;

  while( ip ) {
    if( ip->get_notInList() == true ) {
      if( ip == first ) {
        delete_first();
        ip = first;
      }
      else {
        temp = ip->get_next();
        delete_node( ip );
        ip = temp;
      }
    }
    else   // node stays in list, check next node
      ip = ip->get_next();
  }
}



/*------------------------  restore_notInList()  ---------------------------
|  Sets data member 'notInList' to false for all nodes in list.  Called    |
|  to put edges back in adjacency list after a temporary modification.     |
--------------------------------------------------------------------------*/

void intList::restore_notInList() {
  IntNode * ip;

  for( ip = first; ip; ip = ip->get_next() ) {
    if( ip->get_notInList() )
      ip->set_notInList( false );
  }
}



/*---------------------------  countAdj()  ---------------------------------
|   Counts neighboring vertices that have an id number > x.                |
--------------------------------------------------------------------------*/

int intList::countAdj( int x ) {
  IntNode * ip;
  int count = 0;

  ip = first;
  while( ip != 0 ) {
    if( (ip->get_id_num() > x) && (ip->get_notInList() == false) )
      count++;
    ip = ip->get_next();
  }

  return count;
}



/*-----------------------------  search()  ---------------------------------
|   Returns true if a list contains an intNode with id number x.           |
--------------------------------------------------------------------------*/

bool intList::search( int  x ) {
  IntNode * ip;
  int id;

  for( ip = first; ip; ip = ip->get_next() ) {
    id = ip->get_id_num();
    if( id == x )
      return true;
  }
  
  return false;
}



/*-----------------------------  operator<<  -------------------------------
|   Prints the alleles in the list.                                        |
|   Used for debugging.                                                    |
--------------------------------------------------------------------------*/

ostream & operator<<( ostream & out, const intList & list ) {
  IntNode * currentPtr = list.first;

  if(list.first == 0) {
    out << "The list is empty." << endl;
    return out;
  }

  while(currentPtr != 0) {
    currentPtr->print(out);
    currentPtr = currentPtr->get_next();
  }
  cout << endl;

  return out;
}
