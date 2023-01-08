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
        m_iNumNeighs = stoi( m_vParams[3] );
        arrhenius( stod(m_vParams[ 1 ]), stod(m_vParams[ 2 ]), m_pUtilParams->getTemperature(), (m_iNumNeighs + 1) );
    }
    else if (m_sType.compare("constant") == 0){
        m_dDesorptionRate = stod( m_vParams[1] );
        m_fType = &Desorption::constantType;

        (this->*m_fType)();
    }
    else {
        m_error->error_simple_msg("Not supported type of process -> " + m_sProcName + " | " + m_sType );
        EXIT
    }

    //Create the rule for the adsoprtion process.
    if ( m_bAllNeihs &&  isPartOfGrowth( m_sDesorbed ) )
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

bool Desorption::difSpeciesRule( Site* s){

    //1. Calculate if there are sites at the same height and not oocupied - their number is defined by stoichiometry of the adsorption reaction
    //2. If 1 holds then return true
    //1. Return false

    if ( s->isOccupied() )
        return true;

    return false;
}

void Desorption::constantType(){
    m_dProb = m_dDesorptionRate*m_pLattice->getSize();
}

void Desorption::arrhenius(double v0, double Ed, double T,  int n)
{
    double k = m_pUtilParams->dkBoltz;
    Ed = Ed/m_pUtilParams->dAvogadroNum;

    m_dProb = v0*exp(-(double)n*Ed/(k*T));
}

bool Desorption::rules( Site* s)
{
    (this->*m_fRules)(s);
}

/// This apply for every lattice without a rule.
bool Desorption::allRule( Site* s){
    if ( calculateNeighbors( s ) == m_iNumNeighs )
        return true;
    return false;
}

// This apply for every lattice without a rule which is actually just pick a site and apply it
bool Desorption::basicRule( Site* s){
    return true;
}

void Desorption::perform( Site* s)
{
    (this->*m_fPerform)(s);
}

void Desorption::singleSpeciesSimpleDesorption(Site *s) {
    //For PVD results
    s->decreaseHeight( 1 );
    calculateNeighbors( s ) ;
    m_seAffectedSites.insert( s );
    for ( Site* neigh:s->getNeighs() ) {
        calculateNeighbors( neigh );
        m_seAffectedSites.insert( neigh );

        for ( Site* firstNeigh:neigh->getNeighs() ){
            firstNeigh->setNeighsNum( calculateNeighbors( firstNeigh ) );
            m_seAffectedSites.insert( firstNeigh );
        }
    }
}

void Desorption::multiSpeciesSimpleDesorption(Site *s)
{
    s->setOccupied( false );
    s->setLabel( s->getBelowLabel() );

    m_seAffectedSites.insert( s );
    for ( Site* neigh:s->getNeighs() )
        m_seAffectedSites.insert( neigh );
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
        if ( s->getID() == m_pLattice->getSite( j, m_pLattice->getX() - 1 )->getID() ){
            return true;
        }
    }

    return false;
}

double Desorption::getProbability(){ return m_dProb; }

}
