//============================================================================
//    Apothesis: A kinetic Monte Calro (KMC) code for deposotion processes.
//    Copyright (C) 2019  Nikolaos (Nikos) Cheimarios
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.

//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.

//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <https://www.gnu.org/licenses/>.
//============================================================================
#include "desorption_simple_cubic.h"

namespace MicroProcesses
{

REGISTER_PROCESS_IMPL(Desorption);

Desorption::Desorption():m_iNeigh(0){}
Desorption::~Desorption(){}

bool Desorption::rules( Site* s)
{
    if ( mf_calculateNeighbors( s ) == any_cast<int>(m_mParams["neighs"] ) )
        return true;
    return false;
}

void Desorption::perform( Site* s)
{
    //For PVD results
    s->decreaseHeight( 1 );
    mf_calculateNeighbors( s ) ;
    m_seAffectedSites.insert( s );
    for ( Site* neigh:s->getNeighs() ) {
        mf_calculateNeighbors( neigh );
        m_seAffectedSites.insert( neigh );

        for ( Site* firstNeigh:neigh->getNeighs() ){
            firstNeigh->setNeighsNum( mf_calculateNeighbors( firstNeigh ) );
            m_seAffectedSites.insert( firstNeigh );
        }
    }
}

int Desorption::mf_calculateNeighbors(Site* s)
{

    //We do not need to count the neighbours here!!!
    //We need it only in the rules!
    int neighs = 1;
    for ( Site* neigh:s->getNeighs() ) {
        if ( s->isLowerStep() && neigh->isHigherStep() ){
            if ( neigh->getHeight() >= s->getHeight() + m_pLattice->getStepDiff() + 1 )
                neighs++;
        }
        else if ( neigh->isLowerStep() && s->isHigherStep() ){
            if ( neigh->getHeight() >= s->getHeight() - m_pLattice->getStepDiff() + 1 )
                neighs++;
        }
        else {
            if ( neigh->getHeight() >= s->getHeight() )
                neighs++;
        }
    }

    s->setNeighsNum( neighs );
    return neighs;

    //For flat surfaces
/*    int neighs = 1;
    if ( mf_isInLowerStep( s ) ){
        for ( Site* neigh:s->getNeighs() ) {
            if ( mf_isInHigherStep( neigh ) ){
                if ( neigh->getHeight() >= s->getHeight() + m_pLattice->getStepDiff() + 1 )
                    neighs++;
            }
            else{
                if ( neigh->getHeight() >= s->getHeight() )
                    neighs++;
            }
        }
    }
    else if ( mf_isInHigherStep( s ) ){
        for ( Site* neigh:s->getNeighs() ) {
            if ( mf_isInLowerStep( neigh ) ){
                if ( neigh->getHeight() >= s->getHeight() - (m_pLattice->getStepDiff() + 1 ) )
                    neighs++;
            }
            else{
                if ( neigh->getHeight() >= s->getHeight() )
                    neighs++;
            }
        }
    }
    else {
        for ( Site* neigh:s->getNeighs() ) {
            if ( neigh->getHeight() >= s->getHeight() )
                neighs++;
        }
    }

    s->setNeighsNum( neighs );

    return neighs; */

  //For flat surfaces
/*    int neighs = 1;
    for ( Site* neigh:s->getNeighs() ) {
        if ( neigh->getHeight() >= s->getHeight() )
            neighs++;
    }
    return neighs;*/
}

bool Desorption::mf_isInLowerStep(Site* s)
{
    for (int j = 0; j < m_pLattice->getY(); j++)
        if ( s->getID() == m_pLattice->getSite( j, 0 )->getID() )
            return true;

    return false;
}

bool Desorption::mf_isInHigherStep(Site* s)
{
    for (int j = 0; j < m_pLattice->getY(); j++){
   //     cout<< m_pLattice->getSite( j, m_pLattice->getX() - 1 )->getID() << endl;
        if ( s->getID() == m_pLattice->getSite( j, m_pLattice->getX() - 1 )->getID() ){
            return true;
        }
    }

    return false;
}

double Desorption::getProbability(){

    //These must trenafered in the global definitions
    /*--- Taken from  Lam and Vlachos (2000)PHYSICAL REVIEW B, VOLUME 64, 035401 - DOI: 10.1103/PhysRevB.64.035401 ---*/
    double Na = 6.0221417930e+23;				// Avogadro's number [1/mol]
    double T = any_cast<double>(m_mParams["T"]); //500;						// [K]
    double k = any_cast<double>(m_mParams["k"]); // 1.3806503e-23;			// Boltzmann's constant [j/K]
    double E_d = (7.14e+4)/Na;			// [j]
    double E = 7e+10/Na;  //71128/Na;   	// [j] -> 17 kcal
    double v0 = 1.0e+13;				// [s^-1]
    /*--------------------------------------------------*/

    return v0*exp(-(double)any_cast<int>(m_mParams["neighs"])*E/(k*T));			//DesorptionSimpleCubic 1 neigh
}

}
