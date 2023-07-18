#include "diffusion_rules.h"

bool diffusionBasicRule(Diffusion* proc, Site* s){

    if ( !s->isOccupied() ) return false;

    for ( Site* neigh:s->getNeighs() ) {
        if ( !neigh->isOccupied() && s->getHeight() == neigh->getHeight()  )
            return true;
    }

    return false;
}

bool diffusionAllRule(Diffusion* proc, Site* s){

    if ( proc->calculateNeighbors(s) == proc->getNumNeighs() )
        return true;

    return false;
}
