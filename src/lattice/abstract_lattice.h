#ifndef ABSTRACT_LATTICE_H
#define ABSTRACT_LATTICE_H

#include <string>

using namespace std;
class Lattice;
class Apothesis;

/** The abstract class which is used for the lattice factory **/
class AbstractLattice
  {
  public:
    /// Constructor
    AbstractLattice( const string& );

    /// Destructor
    virtual ~AbstractLattice(){}

    /// Pure virtual method for creating a process.
    virtual Lattice* create(Apothesis* apothesis)= 0; 
  };

#endif // ABSTRACT_LATTICE_H