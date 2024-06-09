#include "abstract_lattice.h"
#include "factory_lattice.h"

AbstractLattice::AbstractLattice(const string& name){
        FactoryLattice::registerThis( name, this);
}