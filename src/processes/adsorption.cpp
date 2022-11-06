//============================================================================
//    Apothesis: A kinetic Monte Calro (KMC) code for deposition processes.
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

Adsorption::Adsorption(){;}
Adsorption::~Adsorption(){;}

void Adsorption::init( vector<string> params )
{
    //Here the params of this process are set and the probability is calcylated (either directly or though calling to a function.
    m_vParams = params;

    //In the first must always be the type```
    m_sType = any_cast<string>(m_vParams[ 0 ]);
    if ( m_sType.compare("simple") == 0 ){
        m_dStick = stod(m_vParams[ 1 ]);
        m_dF = stod(m_vParams[ 2 ]);
        m_dCtot = stod(m_vParams[ 3 ]);
        m_dMW = stod(m_vParams[ 4 ]);

        simple();
    }
    else {
        m_error->error_simple_msg("Not supported type of process: " + m_sType );
        EXIT;
    }
}

void Adsorption::simple()
{
    //These must trenafered in the global definitions
    //      double Na = m_pUtilParams->dAvogadroNum; //6.0221417930e+23;		// Avogadro's number [1/mol]
    //      double P = 1013;					// [Pa]
    //      double T = any_cast<double>(m_vParams[1]); //500;						// [K]
    //      double k = any_cast<double>(m_vParams[2]); // 1.3806503e-23;			// Boltzmann's constant [j/K]
    //      double s0 = any_cast<double>(m_vParams[3]); //0.1;
    //      double C_tot = any_cast<double>(m_vParams[4]);			// [sites/m^2] Vlachos code says [moles sites/m^2]
    //      double m = 32e-3/Na;				// [kg/mol] this is the molecular wei
    //      double y = 2.0e-4;					// Mole fraction of the precursor on the wafer

    double pi = m_pUtilParams->dPi;
    double Na = m_pUtilParams->dAvogadroNum;
    double mass = m_dMW/Na;
    double T = m_pUtilParams->getTemperature();
    double P = m_pUtilParams->getPressure();

    m_dProb = m_dStick*m_dF*P/(m_dCtot*sqrt(2.0e0*pi*mass*m_pUtilParams->dkBoltz*T) );
}


bool Adsorption::rules( Site* s )
{
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

    if (m_pLattice->hasSteps() ){
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
    } else {
        int neighs = 1;
        for ( Site* neigh:s->getNeighs() ) {
            if ( neigh->getHeight() >= s->getHeight() )
                neighs++;
        }
    }

    s->setNeighsNum( neighs );
    return neighs;
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

double Adsorption::getProbability(){ return m_dProb; }

}
