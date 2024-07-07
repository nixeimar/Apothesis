#ifndef REGISTERLATTICE_H
#define REGISTERLATTICE_H

#include "abstract_lattice.h"
#include "lattice.h"

/** A template class used by the factory method for creating the lattice.*/

template<class T>
class Register: public AbstractLattice// Pointers class contains ptrs to master copy of
{
public:
    /// Contructor
    Register<T>(const std::string& name):AbstractLattice( name ){}

    /// Destructor
    virtual ~Register<T>(){}

    /// Creator of the lattice.
    virtual Lattice* create(){
        return new T;
    }

};

#endif // REGISTER_H