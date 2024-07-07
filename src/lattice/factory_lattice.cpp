#include "factory_lattice.h"

FactoryLattice::FactoryLattice(){ }
FactoryLattice::~FactoryLattice(){ }

void FactoryLattice::registerThis( const string& type, AbstractLattice* lattice ){
    getTable()[type] = lattice;
}

Lattice* FactoryLattice::createLattice( const string& type){

    map< string, AbstractLattice* >::iterator it = getTable().find(type);

    if  ( it != getTable().end() )
        return it->second->create();
    else
        return (Lattice*)( 0 );
}

std::map< string, AbstractLattice* >& FactoryLattice::getTable(){
    static std::map<std::string, AbstractLattice*> table;
    return table;
}
