#include "factory_lattice.h"
#include "lattice.h"

using namespace MicroProcesses;

FactoryLattice::FactoryLattice(){ }
FactoryLattice::~FactoryLattice(){ }

void FactoryLattice::registerThis( const string& s, AbstractLattice* proc ){
    getTable()[ s] = proc;
}

Lattice* FactoryLattice::createLattice( const string& s){

    map< string, AbstractLattice* >::iterator it = getTable().find( s);

    if  ( it != getTable().end() )
        return it->second->create();
    else
        return (Lattice*)( 0 );
}

std::map< string, AbstractLattice* >& FactoryLattice::getTable(){
    static std::map<std::string, AbstractLattice*> table;
    return table;
}


