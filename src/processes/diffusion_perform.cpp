#include "diffusion.h"

namespace MicroProcesses
{

void Diffusion::simpleDiffusion( Site* s){

    vector<Site* > toDiffuse;
    for ( Site* neigh:s->getNeighs() ) {
        if ( !neigh->isOccupied() && s->getHeight() == neigh->getHeight() )
            toDiffuse.push_back( neigh );
    }

    // Random pick a site to re-adsorb
    Site* diffuseSite;
    if (  this->getRandomGen() )
        diffuseSite = s->getNeighs().at( this->getRandomGen()->getIntRandom(0, toDiffuse.size()-1 ) );
    else{
        cout << "The random generator has not been defined." << endl;
        EXIT
    }

    diffuseSite->setLabel( s->getLabel() );
    diffuseSite->setOccupied(true);
    s->setLabel( s->getBelowLabel() );
    s->setOccupied(false);

    this->addAffectedSite( diffuseSite ) ;
    for ( Site* neigh:diffuseSite->getNeighs() )
        this->addAffectedSite( neigh ) ;

    this->addAffectedSite( s ) ;
    for ( Site* neigh:s->getNeighs() )
        this->addAffectedSite( neigh );
}


//This is the case of Lam and Vlachos
void Diffusion::performPVD( Site* s){

    //----- This is desorption ------------------------------------------------------------->
    s->decreaseHeight( 1 );
    this->calculateNeighbors( s ) ;
    this->getAffectedSites().insert( s );
    for ( Site* neigh:s->getNeighs() ) {
        this->calculateNeighbors( neigh );
        this->getAffectedSites().insert( neigh );

        for ( Site* firstNeigh:neigh->getNeighs() ){
            firstNeigh->setNeighsNum( this->calculateNeighbors( firstNeigh ) );
            this->getAffectedSites().insert( firstNeigh );
        }
    }
    //--------------------------------------------------------------------------------------<

    // Random pick a site to re-adsorpt
    Site* adsorbSite;
    if (  this->getRandomGen() )
        adsorbSite = s->getNeighs().at( this->getRandomGen()->getIntRandom(0, this->getNumNeighs()-2 ) );
    else{
        cout << "The random generator has not been defined." << endl;
        EXIT
    }

    //----- This is adsoprtion ------------------------------------------------------------->
    s->increaseHeight( 1 );
    this->calculateNeighbors( s );

    this->getAffectedSites().insert( s ) ;

    for ( Site* neigh:s->getNeighs() ) {
        this->calculateNeighbors( neigh );
        this->getAffectedSites().insert( neigh ) ;
    }
    //--------------------------------------------------------------------------------------<
}

}
