#include "diffusion.h"

namespace MicroProcesses
{

bool Diffusion::diffusionBasicRule( Site* s){

    if ( !s->isOccupied() ) return false;

    for ( Site* neigh:s->getNeighs() ) {
        if ( !neigh->isOccupied() && s->getHeight() == neigh->getHeight()  )
            return true;
    }

    return false;
}

bool Diffusion::diffusionAllRule( Site* s){

    if ( calculateNeighbors(s) == this->getNumNeighs() )
        return true;

    return false;
}

}
