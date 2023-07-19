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

#include "desorption.h"

#include "desorption_types.cpp"
#include "desorption_rules.cpp"
#include "desorption_perform.cpp"

namespace MicroProcesses
{

REGISTER_PROCESS_IMPL(Desorption);

Desorption::Desorption():m_bAllNeihs(false){}
Desorption::~Desorption(){}

void Desorption::init(vector<string> params)
{
    //Here the params of this process are set and the probability is calcylated (either directly or though calling to a function.
    m_vParams = params;

    //In the first must always be the type
    m_sType = any_cast<string>(m_vParams[ 0 ]);

    if ( m_sType.compare("arrhenius") == 0 ){
        m_dv0 = stod(m_vParams[ 1 ]);
        m_dEd = stod(m_vParams[ 2 ]);

        m_fType = &Desorption::arrheniusType;
    }
    else if (m_sType.compare("constant") == 0){
        m_dDesorptionRate = stod( m_vParams[1] );

        m_fType = &Desorption::constantType;
    }
    else {
        m_error->error_simple_msg("Not supported type of process -> " + m_sProcName + " | " + m_sType );
        EXIT
    }

    //Set the type of the process
    (this->*m_fType)();

    //Create the rule for the adsoprtion process.
    if ( m_bAllNeihs && isPartOfGrowth( m_sDesorbed ) )
        m_fRules = &Desorption::allRule;
    else if ( !m_bAllNeihs &&  isPartOfGrowth( m_sDesorbed ) )
        m_fRules = &Desorption::basicRule;
    else
        m_fRules = &Desorption::difSpeciesRule;

    //Check what process should be performed.
    //Desorption in PVD will lead to increasing the height of the site
    //Desorption in CVD/ALD will only change the label of the site
    if ( isPartOfGrowth( m_sDesorbed ) )
        m_fPerform = &Desorption::singleSpeciesSimpleDesorption;
    else
        m_fPerform = &Desorption::multiSpeciesSimpleDesorption;
}

bool Desorption::rules( Site* s)
{
    (this->*m_fRules)(s);
}

void Desorption::perform( Site* s)
{
    m_seAffectedSites.clear();
    (this->*m_fPerform)(s);
}

int Desorption::calculateNeighbors(Site* s)
{
    int neighs = 0;

    if ( m_pLattice->hasSteps() ){
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
    }
    else {
        for ( Site* neigh:s->getNeighs() ) {
            if ( neigh->getHeight() >= s->getHeight() )
                neighs++;
        }

        s->setNeighsNum( neighs );
    }

    return neighs;
}

bool Desorption::isInLowerStep(Site* s)
{
    for (int j = 0; j < m_pLattice->getY(); j++)
        if ( s->getID() == m_pLattice->getSite( j, 0 )->getID() )
            return true;

    return false;
}

bool Desorption::isInHigherStep(Site* s)
{
    for (int j = 0; j < m_pLattice->getY(); j++){
        if ( s->getID() == m_pLattice->getSite( j, m_pLattice->getX() - 1 )->getID() )
            return true;
    }

    return false;
}

}
