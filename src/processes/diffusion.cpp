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

#include "diffusion.h"

namespace MicroProcesses
{

REGISTER_PROCESS_IMPL(Diffusion)

Diffusion::Diffusion():m_iNeighNum(0){}
Diffusion::~Diffusion(){}

void Diffusion::perform( int siteID )
{
    m_pLattice->desorp( siteID, m_Species );

    int targetID = m_pLattice->getSite( siteID )->getNeighs().at(  rand()%4  )->getID();

    //From this site get the neighobours id and peformt it.
    m_pLattice->adsorp( targetID, m_Species );
}

double Diffusion::getProbability(){

    //These must trenafered in the global definitions
    /*--- Taken from  Lam and Vlachos (2000)PHYSICAL REVIEW B, VOLUME 64, 035401 - DOI: 10.1103/PhysRevB.64.035401 ---*/
    double Na = 6.0221417930e+23;				// Avogadro's number [1/mol]
    double P = 101325;					// [Pa]
    double T = 500;						// [K]
    double k = 1.3806503e-23;			// Boltzmann's constant [j/K]
    double s0 = 0.1;
    double C_tot = 1.0e+19;				// [sites/m^2] Vlachos code says [moles sites/m^2]
    double E_d = (7.14e+4)/Na;			// [j]
    double E = 71128/Na;   //(7.14e+4)/Na;			// [j] -> 17 kcal
    double m = 32e-3/Na;				// [kg]
    double E_m = (4.28e+4)/Na;			// [j]
    double k_d = 1.0e+13;				// [s^-1]
    double y = 2.0e-3;					// Mole fraction of the precursor on the wafer
    /*--------------------------------------------------*/

    double v0 = k_d; //*exp(-E/(k*T));
    double A = exp( (E_d-E_m)/(k*T) );

    //--------------------- Transitions probability ----------------------------------------//
    return 0; // A*v0*exp( -m_iNeighNum*E/(k*T) );
    //----------------------------------------------------------------------------------------//
}

}
