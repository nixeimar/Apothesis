#include "properties.h"


namespace Utils {

Properties::Properties(Apothesis* apothesis ):
    Pointers( apothesis ),m_dRoughness(0.0), m_dRMS(0.0),m_dEvGrRate(0.0)
{}


double Properties::getMicroroughness()
{
    double dRough = 0.0;
    for ( int i = 0; i < m_lattice->getSize(); i++){
        for (Site* s:m_lattice->getSite( i )->getNeighs() )
            dRough += abs( s->getHeight() - m_lattice->getSite( i )->getHeight() );
    }

    return 1. + dRough/(2.*m_lattice->getSize());
}

double Properties::getRMS()
{
    double sum = 0;
    double mean = 0;
    double dev = 0;

    for (unsigned int i=0; i< m_lattice->getSize(); i++)
        sum += m_lattice->getSite( i )->getHeight();

    mean = sum/(double)m_lattice->getSize();

    for (unsigned int i=0; i<m_lattice->getSize(); i++)
        dev += m_lattice->getSite( i )->getHeight()*m_lattice->getSite( i )->getHeight();

    return sqrt( dev/(double)m_lattice->getSize() );
}

double Properties::eventCountingGrowthRate( int adsorptionCounts, int desortionCounts, double time)
{
    double ra = adsorptionCounts/time;
    double rb = desortionCounts/time;

    return (ra-rb)/(double)m_lattice->getSize(); //growth rate in [ML/s]
}

double Properties::getMeanDH()
{
    double mean = 0.0, sum = 0.0 ;
    if ( m_lattice->getType() == Lattice::FCC ){
        int iCount = 0;
        for (unsigned int i=0; i< m_lattice->getSize(); i++){
            //This is not correct. It should just counts the height. What it is there should be seen by the individual processes.
            if ( m_lattice->getSite( i )->getLabel() == "Cu"){
                if ( m_lattice->getSite(i)->getHeight() > m_lattice->getSite(i)->get1stNeihbors()[ -1 ][ 0 ]->getHeight() ){
                    sum += m_lattice->getSite( i )->getHeight();
                    iCount++;
                }
            }
        }

        mean = sum/iCount;
    }
    else if ( m_lattice->getType() == Lattice::BCC ){
        for (unsigned int i=0; i< m_lattice->getSize(); i++)
            sum += m_lattice->getSite( i )->getHeight();

        mean = sum/m_lattice->getSize();
    }

    return mean;
}

}
