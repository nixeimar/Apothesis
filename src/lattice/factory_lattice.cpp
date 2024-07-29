#include "factory_lattice.h"
#include "lattice.h"

FactoryLattice::FactoryLattice(){ }
FactoryLattice::~FactoryLattice(){ }

void FactoryLattice::registerThis( const string& s, AbstractLattice* proc ){
    getTable()[ s] = proc;
}

Lattice* FactoryLattice::createLattice( const string& s){
    Lattice *lattice = NULL;

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




