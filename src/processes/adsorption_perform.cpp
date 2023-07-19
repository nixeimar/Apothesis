#include "adsorption.h"

namespace MicroProcesses
{

void Adsorption::signleSpeciesAdsorption(Site *s) {
    //Needs check!
    s->increaseHeight( 1 );
    calculateNeighbors( s );
    m_seAffectedSites.insert( s ) ;

    for ( Site* neigh:s->getNeighs() ) {
        calculateNeighbors( neigh );
        m_seAffectedSites.insert( neigh );
    }

    vector<Site*> neighs = s->getNeighs();

    // Because one is already occupied above
    for ( int i = 0 ; i < m_iNumSites-1; i++) {
        int ranNum = m_pRandomGen->getIntRandom( 0,  neighs.size()-1 );
        Site* neigh = neighs[ ranNum ];
        neigh->increaseHeight(1);
        calculateNeighbors( neigh );
        m_seAffectedSites.insert( neigh ) ;

        for ( Site* neigh2:neigh->getNeighs() ) {
            calculateNeighbors( neigh2 );
            m_seAffectedSites.insert( neigh2 );
        }

        neighs.erase( find( neighs.begin(), neighs.end(), neigh ) );
    }
}

void Adsorption::signleSpeciesSimpleAdsorption(Site *s) {
    s->increaseHeight( 1 );
    calculateNeighbors( s );
    m_seAffectedSites.insert( s ) ;

    for ( Site* neigh:s->getNeighs() ) {
        calculateNeighbors( neigh );
        m_seAffectedSites.insert( neigh );
    }
}

void Adsorption::multiSpeciesSimpleAdsorption(Site *s) {
    //Here must hold the previous site in order to appear in case of multiple species forming the growing film
    s->setOccupied( true );
    s->setBelowLabel( s->getLabel() );
    s->setLabel( m_sAdsorbed );

    m_seAffectedSites.insert( s );
    for ( Site* neigh:s->getNeighs() )
        m_seAffectedSites.insert( neigh ) ;
}

void Adsorption::multiSpeciesAdsorption(Site *s) {
    //Here must hold the previous site in order to appear in case of multiple species forming the growing film
    s->setOccupied( true );
    s->setBelowLabel( s->getLabel() );
    s->setLabel( m_sAdsorbed );

    m_seAffectedSites.insert( s );
    for ( Site* neigh:s->getNeighs() )
        m_seAffectedSites.insert( neigh ) ;

    vector<Site*> neighs = s->getNeighs();

    int iNum = 0;
    while (iNum != m_iNumSites-1 ) {
        int ranNum = m_pRandomGen->getIntRandom( 0,  neighs.size()-1 );

        Site* neigh = neighs[ ranNum ];

        if ( !neigh->isOccupied() && neigh->getHeight() == s->getHeight() ) {
            neigh->setOccupied( true );
            neigh->setBelowLabel( neigh->getLabel() );
            neigh->setLabel( m_sAdsorbed );

            m_seAffectedSites.insert( neigh ) ;
            for ( Site* neigh2:neigh->getNeighs() )
                m_seAffectedSites.insert( neigh2 );

            neighs.erase( find( neighs.begin(), neighs.end(), neigh ) );
            iNum++;
        }
        else
            neighs.erase( find( neighs.begin(), neighs.end(), neigh ) );
    }
}

}
