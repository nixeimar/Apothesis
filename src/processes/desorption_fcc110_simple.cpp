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
#include "desorption_fcc110_simple.h"

namespace MicroProcesses
{

REGISTER_PROCESS_IMPL( DesorptionFCC110Simple )

DesorptionFCC110Simple::DesorptionFCC110Simple():m_Species(0), m_iEnableNeighs(4){}

DesorptionFCC110Simple::~DesorptionFCC110Simple(){}

// It can always desorpts if it is adsorbed ...
bool DesorptionFCC110Simple::rules( Site* s )
{
  //  cout << mf_countNeighs( s )  << endl;
  //  cout << any_cast<int>(m_mParams["neighs"] )  << endl;

    if ( s->getHeight() < s->get1stNeihbors()[ -1 ][ 0 ]->getHeight() ||
         s->getHeight() < s->get1stNeihbors()[ -1 ][ 1 ]->getHeight() ||
         s->getHeight() < s->get1stNeihbors()[ -1 ][ 2 ]->getHeight() ||
         s->getHeight() < s->get1stNeihbors()[ -1 ][ 3 ]->getHeight()
         )
        return false;

    if ( mf_countNeighs( s ) == any_cast<int>(m_mParams["neighs"] ) )
        return true;

    return false;
}

void DesorptionFCC110Simple::perform( Site* s )
{
    m_seAffectedSites.clear();

    //Here we assume that the HAMD just desorbed from the site and a Cu atom added in the surface
    s->increaseHeight( 2 );

    for ( int i =0; i < s->get1stNeihbors()[ -1 ].size(); i++)
        m_seAffectedSites.insert( s->get1stNeihbors()[ -1 ][ i ] );

    //In the same level the neighbors must be removed if for some reason there is no couple
    for ( int i =0; i < s->get1stNeihbors()[ 0 ].size(); i++)
        m_seAffectedSites.insert( s->get1stNeihbors()[ 0 ][ i ] );
}

int DesorptionFCC110Simple::mf_countNeighs( Site* s)
{
    int iCount = 4;
    for ( int i =0; i < s->get1stNeihbors()[ 0 ].size(); i++)
        if (  s->get1stNeihbors()[ 0 ][ i ]->getHeight() >= s->getHeight() )
            iCount++;

    return iCount;
}

void DesorptionFCC110Simple::mf_setNeighsNum( Site* s )
{
    int iCount = 4;
    for ( int i =0; i < s->get1stNeihbors()[ 0 ].size(); i++){
        if (  s->get1stNeihbors()[ 0 ][ i ]->getHeight() >= s->getHeight() )
            iCount++;
        s->setNeighsNum( iCount );
    }
}

double DesorptionFCC110Simple::getProbability()
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
