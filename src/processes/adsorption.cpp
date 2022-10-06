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
#include "adsorption.h"

namespace MicroProcesses
{

REGISTER_PROCESS_IMPL( Adsorption )

Adsorption::Adsorption():m_Species(0){}

Adsorption::~Adsorption(){}

bool Adsorption::rules( Site* s )
{
    //You can always adsorb in simple cubic lattices
    return true;
}

void Adsorption::perform( Site* s )
{
    //For PVD results
    s->increaseHeight( 1 );
    mf_calculateNeighbors( s );
    m_seAffectedSites.insert( s ) ;

    for ( Site* neigh:s->getNeighs() ) {
        mf_calculateNeighbors( neigh );
        m_seAffectedSites.insert( neigh ) ;
    }
}

int Adsorption::mf_calculateNeighbors(Site* s)
{
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
    return neighs; // I do not know if we actual need this to be done here ...

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

    return neighs;*/

    //For flat surfaces
    /*    int neighs = 1;
    for ( Site* neigh:s->getNeighs() ) {
        if ( neigh->getHeight() >= s->getHeight() )
            neighs++;
    }
    return neighs;*/
}

bool Adsorption::mf_isInLowerStep(Site* s)
{
    for (int j = 0; j < m_pLattice->getY(); j++)
        if ( s->getID() == m_pLattice->getSite( j, 0 )->getID() )
            return true;

    return false;
}

bool Adsorption::mf_isInHigherStep(Site* s)
{
    for (int j = 0; j < m_pLattice->getY(); j++)
        if ( s->getID() == s->getID() == m_pLattice->getSite( j, m_pLattice->getX() - 1 )->getID() )
            return true;

    return false;
}


double Adsorption::getProbability(){

    //These must trenafered in the global definitions
    double Na = 6.0221417930e+23;		// Avogadro's number [1/mol]
    double P = 101325;					// [Pa]
    double T = any_cast<double>(m_mParams["T"]); //500;						// [K]
    double k = any_cast<double>(m_mParams["k"]); // 1.3806503e-23;			// Boltzmann's constant [j/K]
    double s0 = any_cast<double>(m_mParams["s0"]); //0.1;
    double C_tot = any_cast<double>(m_mParams["C_tot"]);			// [sites/m^2] Vlachos code says [moles sites/m^2]
    // Ctot for copper: 2e-19 (see
    double m = 32e-3/Na;				// [kg/mol] this is the molecular weight
    double y = 2.0e-4;					// Mole fraction of the precursor on the wafer

    return s0*y*P/(C_tot*sqrt(2.0e0*3.14159265*m*k*T) );
}

}
