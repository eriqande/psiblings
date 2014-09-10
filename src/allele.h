// allele.h
// Declaration of the class Allele.

#ifndef ALLELE_H
#define ALLELE_H

#include <iostream>
#include <string.h>
using namespace std;

class Allele
{
  public:
    Allele( string n, int cv, int ct );
    Allele( const Allele & m );
    ~Allele();

    Allele & operator=( const Allele & m );
    bool operator==( const Allele & m ) const;

    void print( ostream & out ) const;
    Allele * copy() const;

    string get_name() const { return name; }
    int get_conversion() const { return conversion; }
    int get_count() const { return count; }
    Allele * get_next() const { return next; }

    void set_count( int ct ) { count = ct; }
    void set_next( Allele * ap ) { next = ap; }

  private:
    string name;
    int conversion;                             // integer alias for allele name
    int count;                                  // # of occurrences of this allele
    Allele * next;
};

#endif
