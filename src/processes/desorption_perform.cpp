#include "desorption.h"

namespace MicroProcesses
{

void Desorption::singleSpeciesSimpleDesorption(Site *s) {
    //For PVD results
    s->decreaseHeight( 1 );
    calculateNeighbors( s ) ;
    m_seAffectedSites.insert( s );
    for ( Site* neigh:s->getNeighs() ) {
        calculateNeighbors( neigh );
        m_seAffectedSites.insert( neigh );

        for ( Site* firstNeigh:neigh->getNeighs() ){
            firstNeigh->setNeighsNum( calculateNeighbors( firstNeigh ) );
            m_seAffectedSites.insert( firstNeigh );
        }
    }
}

void Desorption::multiSpeciesSimpleDesorption(Site *s)
{
    s->setOccupied( false );
    s->setLabel( s->getBelowLabel() );

    m_seAffectedSites.insert( s );
    for ( Site* neigh:s->getNeighs() )
        m_seAffectedSites.insert( neigh );
}

}
