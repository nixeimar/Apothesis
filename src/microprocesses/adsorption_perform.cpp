#include "adsorption_perform.h"

namespace MicroProcesses
{

void signleSpeciesAdsorption(Adsorption* proc, Site *s) {
    //Needs check!
    s->increaseHeight( 1 );
    proc->calculateNeighbors( s );
    proc->addAffectedSite( s ) ;

    for ( Site* neigh:s->getNeighs() ) {
        proc->calculateNeighbors( neigh );
        proc->addAffectedSite( neigh );
    }

    vector<Site*> neighs = s->getNeighs();

    // Because one is already occupied above
    for ( int i = 0 ; i < proc->getNumSites()-1; i++) {
        int ranNum = proc->getRandomGen()->getIntRandom( 0,  neighs.size()-1 );
        Site* neigh = neighs[ ranNum ];
        neigh->increaseHeight(1);
        proc->calculateNeighbors( neigh );
        proc->addAffectedSite( neigh );

        for ( Site* neigh2:neigh->getNeighs() ) {
            proc->calculateNeighbors( neigh2 );
            proc->addAffectedSite( neigh2 );
        }

        neighs.erase( find( neighs.begin(), neighs.end(), neigh ) );
    }
}

void signleSpeciesSimpleAdsorption(Adsorption* proc, Site *s) {
    s->increaseHeight( 1 );
    proc->calculateNeighbors( s );
    proc->addAffectedSite( s ) ;

    for ( Site* neigh:s->getNeighs() ) {
        proc->calculateNeighbors( neigh );
        proc->addAffectedSite( neigh );
    }
}

void multiSpeciesSimpleAdsorption(Adsorption* proc, Site *s) {
    //Here must hold the previous site in order to appear in case of multiple species forming the growing film
    s->setOccupied( true );
    s->setBelowLabel( s->getLabel() );
    s->setLabel( proc->getAdsorbedSpecies() );

    proc->addAffectedSite( s ) ;
    for ( Site* neigh:s->getNeighs() )
        proc->addAffectedSite( neigh ) ;
}

void multiSpeciesAdsorption(Adsorption* proc, Site *s) {
    //Here must hold the previous site in order to appear in case of multiple species forming the growing film
    s->setOccupied( true );
    s->setBelowLabel( s->getLabel() );
    s->setLabel( proc->getAdsorbedSpecies() );

    proc->addAffectedSite( s ) ;
    for ( Site* neigh:s->getNeighs() )
        proc->addAffectedSite( neigh ) ;

    vector<Site*> neighs = s->getNeighs();

    int iNum = 0;
    while (iNum != proc->getNumSites()-1 ) {
        int ranNum = proc->getRandomGen()->getIntRandom( 0,  neighs.size()-1 );

        Site* neigh = neighs[ ranNum ];

        if ( !neigh->isOccupied() && neigh->getHeight() == s->getHeight() ) {
            neigh->setOccupied( true );
            neigh->setBelowLabel( neigh->getLabel() );
            neigh->setLabel( proc->getAdsorbedSpecies() );

            proc->addAffectedSite( neigh ) ;
            for ( Site* neigh2:neigh->getNeighs() )
               proc->addAffectedSite( neigh2 );

            neighs.erase( find( neighs.begin(), neighs.end(), neigh ) );
            iNum++;
        }
        else
            neighs.erase( find( neighs.begin(), neighs.end(), neigh ) );
    }
}

}
