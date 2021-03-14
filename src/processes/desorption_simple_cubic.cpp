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

REGISTER_PROCESS_IMPL(DesorptionSimpleCubic);

DesorptionSimpleCubic::DesorptionSimpleCubic():m_iNeigh(0){}
DesorptionSimpleCubic::~DesorptionSimpleCubic(){}

bool DesorptionSimpleCubic::rules( Site* s)
{
    if ( s->getNeighsNum() == m_iNeigh )
        return true;
    return false;
}

void DesorptionSimpleCubic::perform( Site* s)
{
    //For PVD results
    s->decreaseHeight( 1 );
    s->setNeighsNum( mf_calculateNeighbors( s ) );

    m_seAffectedSites.insert( s );
    for ( Site* neigh:s->getNeighs() ) {
        neigh->setNeighsNum( mf_calculateNeighbors( neigh ) );
        m_seAffectedSites.insert( neigh );

        for ( Site* firstNeigh:neigh->getNeighs() ){
            firstNeigh->setNeighsNum( mf_calculateNeighbors( firstNeigh ) );
            m_seAffectedSites.insert( firstNeigh );
        }
    }
}

int DesorptionSimpleCubic::mf_calculateNeighbors(Site* s)
{
    int neighs = 1;
    for ( Site* neigh:s->getNeighs() ) {
        if ( neigh->getHeight() >= s->getHeight() )
            neighs++;
    }
    return neighs;
}

set<Site*> DesorptionSimpleCubic::getAffectedSites() { return m_seAffectedSites; }


double DesorptionSimpleCubic::getProbability(){

    //These must trenafered in the global definitions
    /*--- Taken from  Lam and Vlachos (2000)PHYSICAL REVIEW B, VOLUME 64, 035401 - DOI: 10.1103/PhysRevB.64.035401 ---*/
    double Na = 6.0221417930e+23;				// Avogadro's number [1/mol]
    double P = 101325;					// [Pa]
    double T = any_cast<double>(m_mParams["T"]); //500;						// [K]
    double k = any_cast<double>(m_mParams["k"]); // 1.3806503e-23;			// Boltzmann's constant [j/K]
    double s0 = 0.1;
    double C_tot = 1.0e+19;				// [sites/m^2] Vlachos code says [moles sites/m^2]
    double E_d = (7.14e+4)/Na;			// [j]
    double E = 71128/Na; //any_cast<double>(m_mParams["E"]); //71128/Na;   //(7.14e+4)/Na;			// [j] -> 17 kcal
    double m = 32e-3/Na;				// [kg]
    double E_m = (4.28e+4)/Na;			// [j]
    double k_d = 1.0e+13;				// [s^-1]
    double y = 2.0e-3;					// Mole fraction of the precursor on the wafera
    /*--------------------------------------------------*/

    m_iNeigh = any_cast<int>(m_mParams["neighs"]);

    double v0 = k_d; //*exp(-E/(k*T));
    double A = 0.0e0; //exp((E_d-E_m)/(k*T));


    return v0*exp(-(double)m_iNeigh*E/(k*T));			//DesorptionSimpleCubic 1 neigh
}

}
