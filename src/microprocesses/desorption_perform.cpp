#include "desorption_perform.h"

namespace MicroProcesses
{

void singleSpeciesSimpleDesorption(Desorption* proc, Site *s) {
    //For PVD results
    s->decreaseHeight( 1 );
    proc->calculateNeighbors( s ) ;
    proc->addAffectedSite( s );
    for ( Site* neigh:s->getNeighs() ) {
        proc->calculateNeighbors( neigh );
        proc->addAffectedSite( neigh );

        for ( Site* firstNeigh:neigh->getNeighs() ){
            firstNeigh->setNeighsNum( proc->calculateNeighbors( firstNeigh ) );
            proc->addAffectedSite( firstNeigh );
        }
    }
}

void multiSpeciesSimpleDesorption(Desorption* proc, Site *s)
{
    s->setOccupied( false );
    s->setLabel( s->getBelowLabel() );

    proc->addAffectedSite( s );
    for ( Site* neigh:s->getNeighs() )
        proc->addAffectedSite( neigh );
}

}
