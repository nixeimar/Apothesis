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
#include "adsorption_fcc1102s_multi.h"

namespace MicroProcesses
{

REGISTER_PROCESS_IMPL( AdsortpionFCC1102SMulti )

AdsortpionFCC1102SMulti::AdsortpionFCC1102SMulti(): m_iEnableNeighs(4)
{
    /* initialize random seed: */
    srand (time(NULL));
}

AdsortpionFCC1102SMulti::~AdsortpionFCC1102SMulti(){}

bool AdsortpionFCC1102SMulti::rules( Site* s )
{
    if ( s->getLabel() != "Cu")
        return false;

    if ( s->getCoupledSite() )
        return true;

    if ( mf_isLowLevelComplete(s) )
    {
        vector<Site* > sites = mf_findPotCouples( s );
        if ( !sites.empty() )
            return true;

/*        if ( !sites.empty() ) {
            Site* coupledSite = sites[ rand() % sites.size() ];
            if ( coupledSite ) {
                s->setCoupledSite( coupledSite );
                coupledSite->setCoupledSite( s );
                return true;
            }
        }*/
    }

    return false;
}

void AdsortpionFCC1102SMulti::perform( Site* s )
{
    m_seAffectedSites.clear();

    //Find its couple
    if ( !s->getCoupledSite() ) {
        vector<Site* > sites = mf_findPotCouples( s );
        Site* coupledSite = sites[ rand() % sites.size() ];
        if ( coupledSite ) {
            s->setCoupledSite( coupledSite );
            coupledSite->setCoupledSite( s );
        }
    }

    s->setLabel("HAMD");
    s->getCoupledSite()->setLabel("HAMD");

    s->increaseHeight( 2 );
    s->getCoupledSite()->increaseHeight( 2 );

    m_seAffectedSites.insert( s );
    m_seAffectedSites.insert( s->getCoupledSite() );

    //Mark the affected sites
    for ( int i =0; i < s->get1stNeihbors()[ -1 ].size(); i++)
        m_seAffectedSites.insert( s->get1stNeihbors()[ -1 ][ i ] );

    for ( int i =0; i < s->get1stNeihbors()[ 0 ].size(); i++)
        m_seAffectedSites.insert( s->get1stNeihbors()[ 0 ][ i ] );

    for ( int i =0; i < s->getCoupledSite()->get1stNeihbors()[ -1 ].size(); i++)
        m_seAffectedSites.insert( s->getCoupledSite()->get1stNeihbors()[ -1 ][ i ] );

    for ( int i =0; i < s->getCoupledSite()->get1stNeihbors()[ 0 ].size(); i++)
        m_seAffectedSites.insert( s->getCoupledSite()->get1stNeihbors()[ 0 ][ i ] );
}

bool AdsortpionFCC1102SMulti::mf_hasCouple( Site* s )
{
    for ( int i = 0; i < s->get1stNeihbors()[ 0 ].size(); i++){
        if (  s->get1stNeihbors()[ 0 ][ i ]->getLabel() == "Cu" && s->getHeight() == s->get1stNeihbors()[ 0 ][ i ]->getHeight()  &&
             mf_isLowLevelComplete(  s->get1stNeihbors()[ 0 ][ i ] ) )
            return true;
    }
    return false;
}

bool AdsortpionFCC1102SMulti::mf_isLowLevelComplete( Site* s )
{
    int iCount = 1;

    //These are the neighbour a level below
    if ( s->get1stNeihbors()[ -1 ][ 0 ]->getHeight() == s->get1stNeihbors()[ -1 ][ 1 ]->getHeight() )
        iCount++;

    if ( s->get1stNeihbors()[ -1 ][ 1 ]->getHeight() == s->get1stNeihbors()[ -1 ][ 2 ]->getHeight() )
        iCount++;

    if ( s->get1stNeihbors()[ -1 ][ 2 ]->getHeight() == s->get1stNeihbors()[ -1 ][ 3 ]->getHeight() )
        iCount++;

    if  ( s->get1stNeihbors()[ -1 ][ 0 ]->getLabel() != "Cu" || s->get1stNeihbors()[ -1 ][ 1 ]->getLabel() != "Cu"
          || s->get1stNeihbors()[ -1 ][ 2 ]->getLabel() != "Cu" || s->get1stNeihbors()[ -1 ][ 3 ]->getLabel() != "Cu")
        return false;

    if ( iCount == m_iEnableNeighs && s->getHeight() < s->get1stNeihbors()[ -1 ][ 0 ]->getHeight() &&  s->getLabel() == "Cu")
        return true;

    return false;
}

vector<Site* > AdsortpionFCC1102SMulti::mf_findPotCouples( Site* s)
{
    vector<Site* > sites;
    for (int i =0; i <  s->get1stNeihbors()[ 0 ].size(); i++)
        //if ( !s->getCoupledSite() && !s->get1stNeihbors()[ 0 ][i]->getCoupledSite() && s->get1stNeihbors()[ 0 ][ i ]->getLabel() == "Cu"
        if ( !s->get1stNeihbors()[ 0 ][i]->getCoupledSite() && s->get1stNeihbors()[ 0 ][ i ]->getLabel() == "Cu"             && mf_isLowLevelComplete( s->get1stNeihbors()[ 0 ][ i ] ) && s->get1stNeihbors()[ 0][ i ]->getHeight() == s->getHeight() )
            sites.push_back( s->get1stNeihbors()[ 0][ i ] );
    return sites;
}

double AdsortpionFCC1102SMulti::getProbability()
{
    //These must trenafered in the global definitions
    double Na = 6.0221417930e+23;		// Avogadro's number [1/mol]
    double P = 1333; //101325;					// [Pa]
    double T = any_cast<double>(m_mParams["T"]); //500;						// [K]
    double k = any_cast<double>(m_mParams["k"]); // 1.3806503e-23;			// Boltzmann's constant [j/K]
    double s0 = any_cast<double>(m_mParams["s0"]); //0.1;
    double C_tot = any_cast<double>(m_mParams["C_tot"]); // [sites/m^2] Vlachos code says [moles sites/m^2]
    double m = 0.4091/Na; //For CuAMD 409 kg/kmol	// [kg/mol] this is the molecular weight
    double y = any_cast<double>(m_mParams["f"]);					// Mole fraction of the precursor on the wafer
    double h = 6.62607004e-34; //m2 kg / s

    double E = 6.3766627287e-19; //j

    return 1.75E-06*(k*T/h)*exp(-40000/(Na*k*T)); // (k*T/h)*s0*y*P/(C_tot*sqrt(2.0e0*3.14159265*m*k*T) )*exp(-40000/(Na*k*T));
}

}
