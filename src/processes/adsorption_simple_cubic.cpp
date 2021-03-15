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
#include "adsorption_simple_cubic.h"

namespace MicroProcesses
{

REGISTER_PROCESS_IMPL( AdsorptionSimpleCubic )

AdsorptionSimpleCubic::AdsorptionSimpleCubic():m_Species(0){}

AdsorptionSimpleCubic::~AdsorptionSimpleCubic(){}

bool AdsorptionSimpleCubic::rules( Site* s )
{
    //You can always adsorb in simple cubic lattices
    return true;
}

void AdsorptionSimpleCubic::perform( Site* s )
{
    //For PVD results
    s->increaseHeight( 1 );
    s->setNeighsNum( mf_calculateNeighbors( s ) );
    m_seAffectedSites.insert( s ) ;

    for ( Site* neigh:s->getNeighs() ) {
        neigh->setNeighsNum( mf_calculateNeighbors( neigh ) );
        m_seAffectedSites.insert( neigh ) ;

        //We do not need this in adsorption
//        for ( Site* firstNeigh:neigh->getNeighs() ){
 //           firstNeigh->setNeighsNum( mf_calculateNeighbors( firstNeigh ) );
  //          m_seAffectedSites.push_back( firstNeigh );
   //     }
    }
}

int AdsorptionSimpleCubic::mf_calculateNeighbors(Site* s)
{
    int neighs = 1;
    for ( Site* neigh:s->getNeighs() ) {
        if ( neigh->getHeight() >= s->getHeight() )
            neighs++;
    }
    return neighs;
}

double AdsorptionSimpleCubic::getProbability(){

    //These must trenafered in the global definitions
    double Na = 6.0221417930e+23;		// Avogadro's number [1/mol]
    double P = 101325;					// [Pa]
    double T = any_cast<double>(m_mParams["T"]); //500;						// [K]
    double k = any_cast<double>(m_mParams["k"]); // 1.3806503e-23;			// Boltzmann's constant [j/K]
    double s0 = any_cast<double>(m_mParams["s0"]); //0.1;
    double C_tot = any_cast<double>(m_mParams["C_tot"]);			// [sites/m^2] Vlachos code says [moles sites/m^2]
    double m = 32e-3/Na;				// [kg/mol] this is the molecular wei
    double y = 2.0e-3;					// Mole fraction of the precursor on the wafer

    return s0*y*P/(C_tot*sqrt(2.0e0*3.14159265*m*k*T) );
}

}
