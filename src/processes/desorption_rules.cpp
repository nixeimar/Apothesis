#include "desorption.h"


namespace MicroProcesses
{

bool Desorption::allRule( Site* s){
    if ( calculateNeighbors( s ) == m_iNumNeighs )
        return true;
    return false;
}

// This apply for every lattice without a rule which is actually just pick a site and apply it
bool Desorption::basicRule( Site* s){
    return true;
}

bool Desorption::difSpeciesRule( Site* s){

    //1. Calculate if there are sites at the same height and not oocupied - their number is defined by stoichiometry of the adsorption reaction
    //2. If 1 holds then return true
    //1. Return false
    if ( s->isOccupied() )
        return true;

    return false;
}

}

