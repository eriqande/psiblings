#include <R.h>
// allele.cc
// Function definitions for class Allele.


#include "allele.h"


/*-----------------------------  constructor  -----------------------------*/

Allele::Allele( string n, int cv, int ct ) {
  name = n;
  conversion = cv;
  count = ct;
  next = 0;
}



/*--------------------------  copy constructor  ---------------------------*/

Allele::Allele( const Allele & a ) {
  name = a.name;
  conversion = a.conversion;
  count = a.count;
  next = a.next;
}



/*------------------------------  destructor  -----------------------------*/

Allele::~Allele() {
  name = "";
  conversion = 0;
  count = 0;
  next = 0;
}



/*----------------------------  assignment  -------------------------------*/

Allele & Allele::operator=( const Allele & a ) {
  if( this == &a ) return *this;

  name = a.name;
  conversion = a.conversion;
  count = a.count;   
  next = a.next;

  return *this;
}



/*-----------------------------  operator==  ------------------------------*/

bool Allele::operator==( const Allele & a ) const {
  if(   ( name == a.name )
     && (conversion == a.conversion)
     && (count == a.count) )
     return true;

  else return false;
}



/*--------------------------------  print()  ------------------------------*/

void Allele::print( ostream & out ) const {
  out << name << ", conv: " << conversion << ", count: " << count << endl;
}



/*--------------------------------  copy()  -------------------------------*/

Allele * Allele::copy() const {
  Allele * ap;
  ap = new Allele( name, conversion, count );

  ap->set_next( next );

  return ap;
}
