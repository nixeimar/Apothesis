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
#include "desorption_fcc110_multi.h"

namespace MicroProcesses
{

REGISTER_PROCESS_IMPL( DesorptionFCC110Multi )

DesorptionFCC110Multi::DesorptionFCC110Multi():m_Species(0), m_iEnableNeighs(4){}

DesorptionFCC110Multi::~DesorptionFCC110Multi(){}

// It can always desorpts if it is adsorbed ...
bool DesorptionFCC110Multi::rules( Site* s )
{
    if ( s->getLabel() == "HAMD" )
        return true;
    return false;
}

void DesorptionFCC110Multi::perform( Site* s )
{
    // Increase the height
//    s->increaseHeight( 2 );
//    s->getCoupledSite()->increaseHeight( 2 );

    //Here we assume that the HAMD just desorbed from the site and a Cu atom added in the surface
    s->setLabel("Cu");
    s->getCoupledSite()->setLabel("Cu");

    m_seAffectedSites.clear();

    for ( int i =0; i < s->get1stNeihbors()[ -1 ].size(); i++)
        m_seAffectedSites.insert( s->get1stNeihbors()[ -1 ][ i ] );

    //In the same level the neighbors must be removed if for some reason there is no couple
    for ( int i =0; i < s->get1stNeihbors()[ 0 ].size(); i++)
        m_seAffectedSites.insert( s->get1stNeihbors()[ 0 ][ i ] );

    for ( int i =0; i < s->getCoupledSite()->get1stNeihbors()[ -1 ].size(); i++)
        m_seAffectedSites.insert( s->getCoupledSite()->get1stNeihbors()[ -1 ][ i ] );

    //In the same level the neighbors must be removed if for some reason there is no couple
    for ( int i =0; i < s->getCoupledSite()->get1stNeihbors()[ 0 ].size(); i++)
        m_seAffectedSites.insert( s->getCoupledSite()->get1stNeihbors()[ 0 ][ i ] );

    //First remove the coupled site
    s->getCoupledSite()->removeCouple();

    //Then the coupled coupled ... else it crashes (as expected!)  ...
    s->removeCouple();
}

int DesorptionFCC110Multi::mf_countNeighs( Site* s)
{
    int iCount = 4;
    for ( int i =0; i < s->get1stNeihbors()[ 0 ].size(); i++)
        if (  s->get1stNeihbors()[ 0 ][ i ]->getLabel() != "HAMD" && s->getHeight() <= s->get1stNeihbors()[ 0 ][ i ]->getHeight()  )
            iCount++;

    return iCount;
}

double DesorptionFCC110Multi::getProbability()
{
    //These must trenafered in the global definitions
    /*--- Taken from  Lam and Vlachos (2000)PHYSICAL REVIEW B, VOLUME 64, 035401 - DOI: 10.1103/PhysRevB.64.035401 ---*/
    double Na = 6.0221417930e+23;				// Avogadro's number [1/mol]
    double T = any_cast<double>(m_mParams["T"]); //500;						// [K]
    double k = any_cast<double>(m_mParams["k"]); // 1.3806503e-23;			// Boltzmann's constant [j/K]
    double E = 90000/Na;  //71128/Na;   //(7.14e+4)/Na;			// [j] -> 17 kcal
    double v0 = 1.0e+13;				// [s^-1]
    /*--------------------------------------------------*/

    return v0*exp(-E/(k*T));			//DesorptionSimpleCubic 1 neigh
}

}
