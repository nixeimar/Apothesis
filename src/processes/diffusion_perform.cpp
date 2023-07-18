#include "diffusion_perform.h"

void simpleDiffusion( Diffusion* proc,  Site* s){

    vector<Site* > toDiffuse;
    for ( Site* neigh:s->getNeighs() ) {
        if ( !neigh->isOccupied() && s->getHeight() == neigh->getHeight() )
            toDiffuse.push_back( neigh );
    }

    // Random pick a site to re-adsorpt
    Site* diffuseSite;
    if (  proc->getRandomGen() )
        diffuseSite = s->getNeighs().at( proc->getRandomGen()->getIntRandom(0, toDiffuse.size()-1 ) );
    else{
        cout << "The random generator has not been defined." << endl;
        EXIT
    }

    diffuseSite->setLabel( s->getLabel() );
    diffuseSite->setOccupied(true);
    s->setLabel( s->getBelowLabel() );
    s->setOccupied(false);

    proc->addAffectedSite( diffuseSite ) ;
    for ( Site* neigh:diffuseSite->getNeighs() )
        proc->addAffectedSite( neigh ) ;

    proc->addAffectedSite( s ) ;
    for ( Site* neigh:s->getNeighs() )
        proc->addAffectedSite( neigh );
}


//This is the case of Lam and Vlachos
void performPVD(Diffusion* proc, Site* s){

    //----- This is desorption ------------------------------------------------------------->
    s->decreaseHeight( 1 );
    proc->calculateNeighbors( s ) ;
    proc->getAffectedSites().insert( s );
    for ( Site* neigh:s->getNeighs() ) {
        proc->calculateNeighbors( neigh );
        proc->getAffectedSites().insert( neigh );

        for ( Site* firstNeigh:neigh->getNeighs() ){
            firstNeigh->setNeighsNum( proc->calculateNeighbors( firstNeigh ) );
            proc->getAffectedSites().insert( firstNeigh );
        }
    }
    //--------------------------------------------------------------------------------------<

    // Random pick a site to re-adsorpt
    Site* adsorbSite;
    if (  proc->getRandomGen() )
        adsorbSite = s->getNeighs().at( proc->getRandomGen()->getIntRandom(0, proc->getNumNeighs()-2 ) );
    else{
        cout << "The random generator has not been defined." << endl;
        EXIT
    }

    //----- This is adsoprtion ------------------------------------------------------------->
    s->increaseHeight( 1 );
    proc->calculateNeighbors( s );

    proc->getAffectedSites().insert( s ) ;

    for ( Site* neigh:s->getNeighs() ) {
        proc->calculateNeighbors( neigh );
        proc->getAffectedSites().insert( neigh ) ;
    }
    //--------------------------------------------------------------------------------------<
}
