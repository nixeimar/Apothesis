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

        m_fType = &Adsorption::mf_simpleType;
    }
    else if (  m_sType.compare("constant") == 0  ) {
        m_fType = &Adsorption::mf_constantType;
        m_dAdsorptionRate = stod(m_vParams[ 1 ]);
    }
    else {
        m_error->error_simple_msg("Not supported type of process: " + m_sType );
        EXIT
    }

    //Create the rule for the adsoprtion process.
    if ( m_iNumSites == 1 && mf_isPartOfGrowth() ){
        setUncoAccepted( true );
        m_fRules = &Adsorption::mf_uncoRule;
    }
    else if ( m_iNumSites > 1 && mf_isPartOfGrowth() ){
        m_fRules = &Adsorption::mf_basicRule;
    }
    else{}

    //Check what process should be performed.
    //Adsorption in PVD will lead to increasing the height of the site
    //Adsorption in CVD will change the label of the site
    if ( mf_isPartOfGrowth() )
        m_fPerform = &Adsorption::mf_performPVD;
    else
        m_fPerform = &Adsorption::mf_performCVDALD;

    (this->*m_fType)();

    cout << m_dProb << endl;
    cout << endl;
}

void Adsorption::mf_constantType(){
    m_dProb = m_dAdsorptionRate*m_pLattice->getSize();
}

bool Adsorption::mf_uncoRule( Site* ){ return true; }

bool Adsorption::mf_basicRule( Site* s){
    if ( mf_calculateNeighbors(s) == m_iNumSites)
        return true;

    return false;
}

//bool Adsorption::mf_difSpeciesRule( Site* s){
    //1. Calculate if there are sites at the same height and not oocupied - their number is defined by stoichiometry of the adsorption reaction
    //2. If 1 holds then return true
    //1. Return false

//}

bool Adsorption::mf_isPartOfGrowth(){
    if (std::find(m_pUtilParams->getGrowthSpecies().begin(), m_pUtilParams->getGrowthSpecies().end(), m_sAdsorbed ) != m_pUtilParams->getGrowthSpecies().end())
        return true;

    return false;
}

void Adsorption::mf_simpleType()
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

//ToDo: To be implemented and checked
void Adsorption::arrheniusType(){;}

//ToDo: To be implemented and checked
void Adsorption::constantType(){; }

bool Adsorption::rules( Site* s )
{
    (this->*m_fRules)(s);
}
void Adsorption::mf_performPVD(Site *s) {
    //For PVD results
    if ( m_iNumSites == 1) {
        s->increaseHeight( 1 );
        mf_calculateNeighbors( s );
        m_seAffectedSites.insert( s ) ;

        for ( Site* neigh:s->getNeighs() ) {
            mf_calculateNeighbors( neigh );
            m_seAffectedSites.insert( neigh );
        }
    }
    else {

        //Needs check!
        s->increaseHeight( 1 );
        mf_calculateNeighbors( s );
        m_seAffectedSites.insert( s ) ;

        for ( Site* neigh:s->getNeighs() ) {
            mf_calculateNeighbors( neigh );
            m_seAffectedSites.insert( neigh );
        }

        vector<Site*> neighs = s->getNeighs();

        for ( int i = 0 ; i < m_iNumSites; i++) {
            int ranNum = m_pRandomGen->getIntRandom( 0,  neighs.size()-1 );
            Site* neigh = s->getNeighs()[ ranNum ];
            neigh->increaseHeight(1);
            mf_calculateNeighbors( neigh );
            m_seAffectedSites.insert( neigh ) ;

            for ( Site* neigh2:neigh->getNeighs() ) {
                mf_calculateNeighbors( neigh2 );
                m_seAffectedSites.insert( neigh2 );
            }

            neighs.erase( find( neighs.begin(), neighs.end(), neighs[ ranNum ] ) );
        }
    }
}

void Adsorption::mf_performCVDALD(Site *s) {
    //Here must hold the previous site in order to appear in case of multiple species forming the growing film
    s->setOccupied( true );
    s->setBelowLabel( s->getLabel() );
    s->setLabel( m_sAdsorbed );

    for ( Site* neigh:s->getNeighs() )
        m_seAffectedSites.insert( neigh ) ;
}

void Adsorption::perform( Site* s )
{
    (this->*m_fPerform)(s);
}

int Adsorption::mf_calculateNeighbors(Site* s)
{
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
        if ( s->getID() == m_pLattice->getSite( j, m_pLattice->getX() - 1 )->getID() )
            return true;

    return false;
}

double Adsorption::getProbability(){ return m_dProb; }

}
