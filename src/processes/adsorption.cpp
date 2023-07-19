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

#include "adsorption_types.cpp"
#include "adsorption_rules.cpp"
#include "adsorption_perform.cpp"

namespace MicroProcesses
{

REGISTER_PROCESS_IMPL( Adsorption )

Adsorption::Adsorption(){}

Adsorption::~Adsorption(){}

void Adsorption::init( vector<string> params )
{
    //Here the params of this process are set and the probability is calculated (either directly or though calling to a function.
    m_vParams = params;

    //Check the type of this adsorption type.
    //Nore: The first argument must always be the type.
    m_sType = any_cast<string>(m_vParams[ 0 ]);
    if ( m_sType.compare("simple") == 0 ){

        m_dStick = stod(m_vParams[ 1 ]);
        m_dF = stod(m_vParams[ 2 ]);
        m_dCtot = stod(m_vParams[ 3 ]);
        m_dMW = stod(m_vParams[ 4 ]);

        m_fType = &Adsorption::simpleType;
    }
    else if (  m_sType.compare("constant") == 0  ) {
        m_dAdsorptionRate = stod(m_vParams[ 1 ]);
        m_fType = &Adsorption::constantType;
    }
    else {
        m_error->error_simple_msg("Not supported type of process: " + m_sType );
        EXIT
    }

    //Assign the type
    (this->*m_fType)();

    //Create the rule for this adsoprtion process.
    if ( m_iNumSites == 1 && isPartOfGrowth( m_sAdsorbed ) ){
        setUncoAccepted( true );
        m_fRules = &Adsorption::uncoRule;
    }
    else if ( m_iNumSites > 1 && isPartOfGrowth( m_sAdsorbed ) )
        m_fRules = &Adsorption::basicRule;
    else if ( m_iNumSites == 1 && !isPartOfGrowth( m_sAdsorbed ) )
        m_fRules = &Adsorption::multiSpeciesSimpleRule;
    else if ( m_iNumSites > 1 && !isPartOfGrowth( m_sAdsorbed ) )
        m_fRules = &Adsorption::multiSpeciesRule;
    else {
        m_error->error_simple_msg("The rule for this process has not been defined.");
        EXIT
    }

    //Check what process should be performed.
    //Adsorption in PVD will lead to increasing the height of the site
    //Adsorption in CVD/ALD will only change the label of the site. The height will change from surface reaction.
    if ( m_iNumSites == 1  && isPartOfGrowth(m_sAdsorbed) )
        m_fPerform = &Adsorption::signleSpeciesSimpleAdsorption;
    else if ( m_iNumSites > 1  && isPartOfGrowth( m_sAdsorbed ) )
        m_fPerform = &Adsorption::signleSpeciesAdsorption;
    else if ( m_iNumSites == 1 && !isPartOfGrowth(m_sAdsorbed) )
        m_fPerform = &Adsorption::multiSpeciesSimpleAdsorption;
    else if ( m_iNumSites > 1 && !isPartOfGrowth(m_sAdsorbed) )
        m_fPerform = &Adsorption::multiSpeciesAdsorption;
    else {
        m_error->error_simple_msg("The process is not defined | " + m_sProcName );
        EXIT
    }
}

int Adsorption::countVacantSites( Site* s){
    int iCount = 0;
    for (Site* neigh:s->getNeighs() ){
        if ( !neigh->isOccupied() && s->getHeight() == neigh->getHeight() )
            iCount++;
    }

    return iCount;
}

bool Adsorption::rules( Site* s )
{
    (this->*m_fRules)(s);
}

void Adsorption::perform( Site* s ) {

    m_seAffectedSites.clear();
    (this->*m_fPerform)(s);
}

int Adsorption::calculateNeighbors(Site* s){

    int neighs = 0;

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
        int neighs = 0;
        for ( Site* neigh:s->getNeighs() ) {
            if ( neigh->getHeight() >= s->getHeight() )
                neighs++;
        }
    }

    s->setNeighsNum( neighs );
    return neighs;
}

bool Adsorption::isInLowerStep(Site* s) {

    for (int j = 0; j < m_pLattice->getY(); j++)
        if ( s->getID() == m_pLattice->getSite( j, 0 )->getID() )
            return true;

    return false;
}

bool Adsorption::isInHigherStep(Site* s){

    for (int j = 0; j < m_pLattice->getY(); j++)
        if ( s->getID() == m_pLattice->getSite( j, m_pLattice->getX() - 1 )->getID() )
            return true;

    return false;
}

}
