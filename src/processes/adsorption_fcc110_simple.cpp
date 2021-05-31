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
#include "adsorption_fcc110_simple.h"

namespace MicroProcesses
{

REGISTER_PROCESS_IMPL( AdsorptionFCC110Simple )

AdsorptionFCC110Simple::AdsorptionFCC110Simple():m_Species(0){}

AdsorptionFCC110Simple::~AdsorptionFCC110Simple(){}

bool AdsorptionFCC110Simple::rules( Site* s )
{
    if ( s == targetSite )
        return false;

    int iCount = 1;

    //These are the neighbour a level below (which are the same as the level above!)
    if ( s->get1stNeihbors()[ -1 ][ 0 ]->getHeight() == s->get1stNeihbors()[ -1 ][ 1 ]->getHeight() )
        iCount++;

    if ( s->get1stNeihbors()[ -1 ][ 1 ]->getHeight() == s->get1stNeihbors()[ -1 ][ 2 ]->getHeight() )
        iCount++;

    if ( s->get1stNeihbors()[ -1 ][ 2 ]->getHeight() == s->get1stNeihbors()[ -1 ][ 3 ]->getHeight() )
        iCount++;

    if ( iCount == m_iEnableNeighs )
        return true;

    return false;
}

void AdsorptionFCC110Simple::perform( Site* s )
{
    //Set the target site to use it later in the rules
    targetSite = s;
    s->increaseHeight( 2 ); //Since this is FCC

    //Get the neighs in the same level
    for (Site* s:s->get1stNeihbors().at(0) )
        m_seAffectedSites.insert( s );

    //Loop over the sites in the same level
    //The 4 neighs below are already occupied so we start from 4
    int iCount = 4;
    for (Site* s1:s->get1stNeihbors()[ 0 ])
        if ( s1->getHeight() >= s->getHeight() )
            iCount++;
    s->setNeighsNum( iCount );

    //These are the neighbour a level below (which are the same as the level above!)
    for (Site* s1:s->get1stNeihbors()[ -1 ])
        m_seAffectedSites.insert( s1 );
}

double AdsorptionFCC110Simple::getProbability()
{

    //These must trenafered in the global definitions
    double Na = 6.0221417930e+23;		// Avogadro's number [1/mol]
    double P = 101325;					// [Pa]
    double T = any_cast<double>(m_mParams["T"]); //500;						// [K]
    double k = any_cast<double>(m_mParams["k"]); // 1.3806503e-23;			// Boltzmann's constant [j/K]
    double s0 = any_cast<double>(m_mParams["s0"]); //0.1;
    double C_tot = any_cast<double>(m_mParams["C_tot"]);			// [sites/m^2] Vlachos code says [moles sites/m^2]
    double m = 32e-3/Na;				// [kg/mol] this is the molecular wei
    double y = 2.0e-4;					// Mole fraction of the precursor on the wafer

    return s0*y*P/(C_tot*sqrt(2.0e0*3.14159265*m*k*T) );
}

}
