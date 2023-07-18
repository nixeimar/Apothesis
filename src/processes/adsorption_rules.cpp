#include "adsorption_rules.h"

bool basicRule(Adsorption* proc, Site* s){

 //   if ( calculateNeighbors(s) == m_iNumSites)
    if ( proc->calculateNeighbors( s ) == proc->getNumSites() )
        return true;

    return false;
}

bool multiSpeciesSimpleRule(Adsorption* proc, Site* s){
    //1. If the species is not occupied return true
    //2. Return false
    if ( !s->isOccupied() )
        return true;

    return false;
}

bool multiSpeciesRule(Adsorption* proc, Site* s){

    //1. If the species is not occupied
    //2. and if neighbours equal to the sites needed by m_iNumSites are vacant
    //3. and have the same height return true
    //4. Return false
    if ( s->isOccupied() || proc->countVacantSites( s ) != proc->getNumVacantSites() )
        return false;

    return true;
}


