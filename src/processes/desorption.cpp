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
        m_fType = &Desorption::mf_constantType;
    }
    else {
        m_error->error_simple_msg("Not supported type of process -> " + m_sProcName + " | " + m_sType );
        EXIT
    }

    //Create the rule for the adsoprtion process.
    if ( m_bAllNeihs )
        m_fRules = &Desorption::mf_allRule;
    else
        m_fRules = &Desorption::mf_basicRule;

    //Check what process should be performed.
    //Desorption in PVD will lead to increasing the height of the site
    //Desorption in CVD will change the label of the site
    if ( mf_isPartOfGrowth() )
        m_fPerform = &Desorption::mf_performPVD;
    else
        m_fPerform = &Desorption::mf_performCVDALD;

}

bool Desorption::mf_isPartOfGrowth(){
    if (std::find(m_pUtilParams->getGrowthSpecies().begin(), m_pUtilParams->getGrowthSpecies().end(), m_sDesorbed ) != m_pUtilParams->getGrowthSpecies().end())
        return true;

    return false;
}

void Desorption::mf_constantType(){
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
bool Desorption::mf_allRule( Site* s){
    if ( mf_calculateNeighbors( s ) == m_iNumNeighs )
        return true;
    return false;
}

// This apply for every lattice without a rule.
bool Desorption::mf_basicRule( Site* s){
    return true;
}

void Desorption::perform( Site* s)
{
    (this->*m_fPerform)(s);
}

void Desorption::mf_performPVD(Site *s) {
    //For PVD results
    s->decreaseHeight( 1 );
    mf_calculateNeighbors( s ) ;
    m_seAffectedSites.insert( s );
    for ( Site* neigh:s->getNeighs() ) {
        mf_calculateNeighbors( neigh );
        m_seAffectedSites.insert( neigh );

        for ( Site* firstNeigh:neigh->getNeighs() ){
            firstNeigh->setNeighsNum( mf_calculateNeighbors( firstNeigh ) );
            m_seAffectedSites.insert( firstNeigh );
        }
    }
}

void Desorption::mf_performCVDALD(Site *s)
{
    s->setOccupied( false );
    s->setLabel( s->getBelowLabel() );

    m_seAffectedSites.insert( s );
    for ( Site* neigh:s->getNeighs() )
        m_seAffectedSites.insert( neigh );
}

int Desorption::mf_calculateNeighbors(Site* s)
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

bool Desorption::mf_isInLowerStep(Site* s)
{
    for (int j = 0; j < m_pLattice->getY(); j++)
        if ( s->getID() == m_pLattice->getSite( j, 0 )->getID() )
            return true;

    return false;
}

bool Desorption::mf_isInHigherStep(Site* s)
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
