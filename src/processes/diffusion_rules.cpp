#include "diffusion_rules.h"

namespace MicroProcesses
{

bool diffusionBasicRule( Diffusion* proc, Site* s){

    if ( !s->isOccupied() || s->getLabel().compare( proc->getDiffused() ) != 0 ) return false;

    for ( Site* neigh:s->getNeighs() ) {
        if ( !neigh->isOccupied() && s->getHeight() == neigh->getHeight() )
            return true;
    }

    return false;
}

bool diffusionBasicAllRule( Diffusion* proc, Site* s){

    if ( !s->isOccupied() || s->getLabel().compare( proc->getDiffused() ) != 0
         || proc->calculateNeighbors(s) != proc->getNumNeighs() ) return false;

    for ( Site* neigh:s->getNeighs() ) {
        if ( !neigh->isOccupied() && s->getHeight() == neigh->getHeight() )
            return true;
    }

    return false;
}


bool diffusionAllRule( Diffusion* proc, Site* s){

    if ( proc->calculateNeighbors(s) == proc->getNumNeighs() )
        return true;

    return false;
}

}
