#include "desorption_rules.h"

namespace MicroProcesses
{

bool allRule(Desorption* proc, Site* s){
    if ( proc->calculateNeighbors( s ) == proc->getNumNeighs() )
        return true;
    return false;
}

// This apply for every lattice without a rule which is actually just pick a site and apply it
bool basicRule(Desorption* proc, Site* s){
    return true;
}

bool difSpeciesRule(Desorption* proc, Site* s){

    //1. Calculate if there are sitess at the same height and not oocupied - their number is defined by stoichiometry of the adsorption reaction
    //2. If 1 holds then return true
    //1. Return false
    if ( s->isOccupied() )
        return true;

    return false;
}

}

