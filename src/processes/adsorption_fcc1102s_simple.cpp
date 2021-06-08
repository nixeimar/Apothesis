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
#include "adsorption_fcc1102s_simple.h"

namespace MicroProcesses
{

REGISTER_PROCESS_IMPL( AdsorptionFCC1102SSimple )

AdsorptionFCC1102SSimple::AdsorptionFCC1102SSimple(): m_iEnableNeighs(4){}

AdsorptionFCC1102SSimple::~AdsorptionFCC1102SSimple(){}

bool AdsorptionFCC1102SSimple::rules( Site* s )
{
    //This is needed in order to remove the site from the process
    if ( (m_psTargetSite && s == m_psTargetSite ) || (m_psCoupledSite && s == m_psCoupledSite )  )
       return false;

    if ( !m_lToRemove.empty() && std::find(m_lToRemove.begin(), m_lToRemove.end(), s) != m_lToRemove.end() )
        return false;

    int iCountA = 1;
    bool firstCond = false;
    //These are the neighbour a level below (which are the same as the level above!)
    if ( s->get1stNeihbors()[ -1 ][ 0 ]->getHeight() == s->get1stNeihbors()[ -1 ][ 1 ]->getHeight() )
        iCountA++;

    if ( s->get1stNeihbors()[ -1 ][ 1 ]->getHeight() == s->get1stNeihbors()[ -1 ][ 2 ]->getHeight() )
        iCountA++;

    if ( s->get1stNeihbors()[ -1 ][ 2 ]->getHeight() == s->get1stNeihbors()[ -1 ][ 3 ]->getHeight() )
        iCountA++;

    if ( iCountA == m_iEnableNeighs && s->getHeight() <  s->get1stNeihbors()[ -1 ][ 0 ]->getHeight() && mf_hasCouple( s ) )
        return true;

    return false;
}

void AdsorptionFCC1102SSimple::perform( Site* s )
{
    m_seAffectedSites.clear();
    m_lToRemove.clear();

    //Set the target site to use it later in the rules
    m_psTargetSite = s;
    m_seAffectedSites.insert( m_psTargetSite );

    m_psCoupledSite = mf_findCouple();
    m_seAffectedSites.insert( m_psCoupledSite );

    m_psTargetSite->increaseHeight( 2 );
    m_psCoupledSite->increaseHeight( 2);

    //These are the neighbour a level below (which are the same as the level above!) and must be checked
    for ( int i =0; i < m_psTargetSite->get1stNeihbors()[ -1 ].size(); i++)
        m_seAffectedSites.insert( m_psTargetSite->get1stNeihbors()[ -1 ][ i ] );

    for ( int i =0; i < m_psCoupledSite->get1stNeihbors()[ -1 ].size(); i++)
        m_seAffectedSites.insert( m_psCoupledSite->get1stNeihbors()[ -1 ][ i ] );

    //In the same level the neighbors must be removed if for some reason there is no couple
    for ( int i =0; i < m_psTargetSite->get1stNeihbors()[ 0 ].size(); i++){
        if ( !mf_hasCouple( m_psTargetSite->get1stNeihbors()[ 0 ][ i ] ) ) {
            m_lToRemove.push_back( m_psTargetSite->get1stNeihbors()[ 0 ][ i ] );
            m_seAffectedSites.insert( m_psTargetSite->get1stNeihbors()[ 0 ][ i ] );
        }
    }

    for ( int i =0; i < m_psCoupledSite->get1stNeihbors()[ 0 ].size(); i++){
        if ( !mf_hasCouple( m_psCoupledSite->get1stNeihbors()[ 0 ][ i ] ) ) {
            m_lToRemove.push_back( m_psCoupledSite->get1stNeihbors()[ 0 ][ i ] );
            m_seAffectedSites.insert( m_psCoupledSite->get1stNeihbors()[ 0 ][ i ] );
        }
    }
}

bool AdsorptionFCC1102SSimple::mf_hasCouple( Site* s )
{
    for ( int i =0; i < s->get1stNeihbors()[ 0 ].size(); i++){
        if ( s->getHeight() == s->get1stNeihbors()[ 0 ][ i ]->getHeight()  &&
             mf_isLowLevelComplete(  s->get1stNeihbors()[ 0 ][ i ] ) )
            return true;
    }
    return false;
}

bool AdsorptionFCC1102SSimple::mf_isLowLevelComplete( Site* s )
{
    int iCount = 1;
    //These are the neighbour a level below (which are the same as the level above!)
    if ( s->get1stNeihbors()[ -1 ][ 0 ]->getHeight() == s->get1stNeihbors()[ -1 ][ 1 ]->getHeight() )
        iCount++;

    if ( s->get1stNeihbors()[ -1 ][ 1 ]->getHeight() == s->get1stNeihbors()[ -1 ][ 2 ]->getHeight() )
        iCount++;

    if ( s->get1stNeihbors()[ -1 ][ 2 ]->getHeight() == s->get1stNeihbors()[ -1 ][ 3 ]->getHeight() )
        iCount++;

    if ( iCount == m_iEnableNeighs && s->getHeight() <  s->get1stNeihbors()[ -1 ][ 0 ]->getHeight() )
        return true;

    return false;
}


Site* AdsorptionFCC1102SSimple::mf_findCouple()
{
    for (int i =0; i <  m_psTargetSite->get1stNeihbors()[ 0 ].size(); i++)
        if ( m_psTargetSite->get1stNeihbors()[ 0 ][ i ]->getHeight() == m_psTargetSite->getHeight() && mf_isLowLevelComplete( m_psTargetSite->get1stNeihbors()[ 0 ][ i ] ) )
            return m_psTargetSite->get1stNeihbors()[ 0 ][ i ];

    return 0;
}

// Reconsider this ... Can we do it faster e.g. do not count the neighs everytime ????
void AdsorptionFCC1102SSimple::mf_setNeighsNum( Site* s)
{
    int iCount = 4;
    for ( int i =0; i < s->get1stNeihbors()[ 0 ].size(); i++){
        if (  s->get1stNeihbors()[ 0 ][ i ]->getHeight() >= s->getHeight() )
            iCount++;
        s->setNeighsNum( iCount );
    }
}

double AdsorptionFCC1102SSimple::getProbability()
{

    //These must trenafered in the global definitions
    double Na = 6.0221417930e+23;		// Avogadro's number [1/mol]
    double P = 101325;					// [Pa]
    double T = any_cast<double>(m_mParams["T"]); //500;						// [K]
    double k = any_cast<double>(m_mParams["k"]); // 1.3806503e-23;			// Boltzmann's constant [j/K]
    double s0 = any_cast<double>(m_mParams["s0"]); //0.1;
    double C_tot = any_cast<double>(m_mParams["C_tot"]);			// [sites/m^2] Vlachos code says [moles sites/m^2]
    double m = 32e-3/Na;				// [kg/mol] this is the molecular wei
    double y = 2.0e-4;					// Mole fraction of the precursor on the wafer

    return s0*y*P/(C_tot*sqrt(2.0e0*3.14159265*m*k*T) );
}

}
