#ifndef FACTORY_LATTICE_H
#define FACTORY_LATTICE_H

#include <string>

#include "register.h"

class Lattice;
class AbstractLattice;
class Apothesis;

using namespace std;

class FactoryLattice
{
    public:
    /// Constructor
    FactoryLattice();

    /// Destructor
    virtual ~FactoryLattice();

    /// Register the new lattice (Store it in the map table since a map table cannot be initialized)
    static void registerThis( const string& , AbstractLattice* );

    /// Create the lattice
    static Lattice* createLattice( const string&, Apothesis* );

//should we include the getTable for the lattice?

}

#define REGISTER_LATTICE( __NAME__ ) \
    private: static const Register< __NAME__> creator;

#define REGISTER_LATTICE_IMPL( __NAME__) \
    const Register< __NAME__ > __NAME__::creator(#__NAME__);

#endif // FACTORY_LATTICE_H