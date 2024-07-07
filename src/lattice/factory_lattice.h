#ifndef FACTORY_LATTICE_H
#define FACTORY_LATTICE_H

#include <string>
#include <map>
#include "abstract_lattice.h"
#include "register_lattice.h"

class FactoryLattice
{
public:
    /// Constructor
    FactoryLattice();

    /// Destructor
    virtual ~FactoryLattice();

    /// Register the new process (Store it in the map table since a map table cannot be initialized)
    static void registerThis( const string& , AbstractLattice* );

    /// Create the process
    static Lattice* createLattice( const string& );

private:
    /// The factory map
    static map< string, AbstractLattice* >& getTable();

};

#define REGISTER_LATTICE( __NAME__ ) \
    private: static const Register< __NAME__> creator;

#define REGISTER_LATTICE_IMPL( __NAME__) \
    const Register< __NAME__ > __NAME__::creator(#__NAME__);

#endif // FACTORY_LATTICE_H
