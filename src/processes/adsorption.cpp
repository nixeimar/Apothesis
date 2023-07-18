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

Adsorption::Adsorption(){}

Adsorption::~Adsorption(){}

void Adsorption::init( vector<string> params )
{
    //Here the params of this process are set and the probability is calcylated (either directly or though calling to a function.
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

    //Create the rule for the adsoprtion process.
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
    //Adsorption in CVD/ALD will only change the label of the site
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

    (this->*m_fType)();
}

void Adsorption::constantType(){
    m_dProb = m_dAdsorptionRate*m_iNumVacant;
}

bool Adsorption::uncoRule( Site* ){ return true; }

bool Adsorption::basicRule( Site* s){
    if ( calculateNeighbors(s) == m_iNumSites)
        return true;

    return false;
}

bool Adsorption::multiSpeciesSimpleRule( Site* s){
    //1. If the species is not occupied return true
    //2. Return false
    if ( !s->isOccupied() )
        return true;

    return false;
}

bool Adsorption::multiSpeciesRule( Site* s){

    //1. If the species is not occupied
    //2. and if neighbours equal to the sites needed by m_iNumSites are vacant
    //3. and have the same height return true
    //4. Return false
    if ( s->isOccupied() || countVacantSites(s) != m_iNumVacant )
        return false;

    return true;
}

int Adsorption::countVacantSites( Site* s){
    int iCount = 0;
    for (Site* neigh:s->getNeighs() ){
        if ( !neigh->isOccupied() && s->getHeight() == neigh->getHeight() )
            iCount++;
    }

    return iCount;
}

void Adsorption::simpleType()
{
    double pi = m_pUtilParams->dPi;
    double Na = m_pUtilParams->dAvogadroNum; // Avogadro's number [1/mol]
    double mass = m_dMW/Na; //[kg/mol]
    double T = m_pUtilParams->getTemperature(); //[K]
    double P = m_pUtilParams->getPressure(); //[Pa]

    m_dProb = m_dStick*m_dF*P/(m_dCtot*sqrt(2.0e0*pi*mass*m_pUtilParams->dkBoltz*T) );
}

//ToDo: To be implemented and checked
void Adsorption::arrheniusType(){}

bool Adsorption::rules( Site* s )
{
    (this->*m_fRules)(s);
}

void Adsorption::signleSpeciesAdsorption(Site *s) {

    //Needs check!
    s->increaseHeight( 1 );
    calculateNeighbors( s );
    m_seAffectedSites.insert( s ) ;

    for ( Site* neigh:s->getNeighs() ) {
        calculateNeighbors( neigh );
        m_seAffectedSites.insert( neigh );
    }

    vector<Site*> neighs = s->getNeighs();

    // Because one is already occupied above
    for ( int i = 0 ; i < m_iNumSites-1; i++) {
        int ranNum = m_pRandomGen->getIntRandom( 0,  neighs.size()-1 );
        Site* neigh = neighs[ ranNum ];
        neigh->increaseHeight(1);
        calculateNeighbors( neigh );
        m_seAffectedSites.insert( neigh ) ;

        for ( Site* neigh2:neigh->getNeighs() ) {
            calculateNeighbors( neigh2 );
            m_seAffectedSites.insert( neigh2 );
        }

        neighs.erase( find( neighs.begin(), neighs.end(), neigh ) );
    }

    cout << endl;
}

void Adsorption::signleSpeciesSimpleAdsorption(Site *s) {

    s->increaseHeight( 1 );
    calculateNeighbors( s );
    m_seAffectedSites.insert( s ) ;

    for ( Site* neigh:s->getNeighs() ) {
        calculateNeighbors( neigh );
        m_seAffectedSites.insert( neigh );
    }
}

void Adsorption::multiSpeciesSimpleAdsorption(Site *s) {

    //Here must hold the previous site in order to appear in case of multiple species forming the growing film
    s->setOccupied( true );
    s->setBelowLabel( s->getLabel() );
    s->setLabel( m_sAdsorbed );

    m_seAffectedSites.insert( s );
    for ( Site* neigh:s->getNeighs() )
        m_seAffectedSites.insert( neigh ) ;

    cout << endl;
}

void Adsorption::multiSpeciesAdsorption(Site *s) {

    //Here must hold the previous site in order to appear in case of multiple species forming the growing film
    s->setOccupied( true );
    s->setBelowLabel( s->getLabel() );
    s->setLabel( m_sAdsorbed );

    m_seAffectedSites.insert( s );
    for ( Site* neigh:s->getNeighs() )
        m_seAffectedSites.insert( neigh ) ;

    vector<Site*> neighs = s->getNeighs();

    int iNum = 0;
    while (iNum != m_iNumSites-1 ) {
        int ranNum = m_pRandomGen->getIntRandom( 0,  neighs.size()-1 );

        Site* neigh = neighs[ ranNum ];

        if ( !neigh->isOccupied() && neigh->getHeight() == s->getHeight() ) {
            neigh->setOccupied( true );
            neigh->setBelowLabel( neigh->getLabel() );
            neigh->setLabel( m_sAdsorbed );

            m_seAffectedSites.insert( neigh ) ;
            for ( Site* neigh2:neigh->getNeighs() )
                m_seAffectedSites.insert( neigh2 );

            neighs.erase( find( neighs.begin(), neighs.end(), neigh ) );
            iNum++;
        }
        else
            neighs.erase( find( neighs.begin(), neighs.end(), neigh ) );
    }
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

double Adsorption::getRateConstant(){ return m_dProb; }

}
