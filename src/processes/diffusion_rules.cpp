#include "diffusion_rules.h"

namespace MicroProcesses
{

bool diffusionBasicRule( Diffusion* proc, Site* s){

    if ( !s->isOccupied() || s->getLabel().compare( proc->getDiffused() ) != 0 ) return false;

    for ( Site* neigh:s->getNeighs() )
        if ( !neigh->isOccupied() && neigh->getHeight() == s->getHeight() )
            return true;

    return false;
}

bool diffusionBasicAllRule( Diffusion* proc, Site* s){

    if ( !s->isOccupied() || s->getLabel().compare( proc->getDiffused() ) != 0
         || proc->countVacantSites(s) != proc->getNumVacantSites() ) return false;

    for ( Site* neigh:s->getNeighs() )
        if ( !neigh->isOccupied() && neigh->getHeight() == s->getHeight() )
            return true;

    return false;
}


bool diffusionAllRule( Diffusion* proc, Site* s){

    if ( s->isOccupied() || !proc->isPartOfGrowth( s->getLabel() ) ||
         proc->countVacantSites(s) != proc->getNumVacantSites() ) return false;

    for ( Site* neigh:s->getNeighs() )
        if ( !neigh->isOccupied() && s->getHeight() == neigh->getHeight() - 1 )
            return true;

    return false;
}

}
