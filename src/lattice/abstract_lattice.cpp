#include "abstract_lattice.h"
#include "factory_lattice.h"

AbstractLattice::AbstractLattice( const string& type ){
    FactoryLattice::registerThis( type, this);
}