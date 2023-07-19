#include "adsorption.h"

namespace MicroProcesses
{

bool Adsorption::basicRule(Site* s){

    if ( calculateNeighbors(s) == m_iNumSites)
        return true;

    return false;
}

bool Adsorption::multiSpeciesSimpleRule( Site* s){
    //1. If the species is not occupied return true
    //2. Return false
    if ( !s->isOccupied() )
        return true;

    return false;
}

bool Adsorption::multiSpeciesRule( Site* s){

    //1. If the species is not occupied
    //2. and if neighbours equal to the sites needed by m_iNumSites are vacant
    //3. and have the same height return true
    //4. Return false
    if ( s->isOccupied() || countVacantSites(s) != m_iNumVacant )
        return false;

    return true;
}

}


