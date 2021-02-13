#include "properties.h"


namespace Utils {

Properties::Properties(Apothesis* apothesis ):
    Pointers( apothesis ),m_dRoughness(0.0), m_dRMS(0.0),m_dEvGrRate(0.0)
{;}

void Properties::calculateRoughness()
{
    double dRough = 0.0;
    for ( int i = 0; i < m_lattice->getSize(); i++){
        for (Site* s:m_lattice->getSite( i )->getNeighs() )
            dRough += abs( s->getHeight() - m_lattice->getSite( i )->getHeight() );
    }

    m_dRoughness = 1. + dRough/(2.*m_lattice->getSize());
}


void Properties::RMS()
{
    ;
}

void Properties::eventCountingGrowthRate()
{
    ;
}

}
