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


#include "adsroption_simple_cubic_4s_multi.h"

namespace MicroProcesses
{

REGISTER_PROCESS_IMPL( AdsroptionSimpleCubic4sMulti )

AdsroptionSimpleCubic4sMulti::AdsroptionSimpleCubic4sMulti()
{
    /* initialize random seed: */
    //srand(2000);
    srand (time(NULL));
}

AdsroptionSimpleCubic4sMulti::~AdsroptionSimpleCubic4sMulti(){}

bool AdsroptionSimpleCubic4sMulti::rules( Site* s )
{
    if ( s->getLabel() != "Cu") return false;
    m_site = s;

    if ( mf_westUpTest() )
        return true;

    if ( mf_westDownTest() )
        return true;

    if ( mf_eastUpTest() )
        return true;

    if ( mf_eastDownTest() )
        return true;

    return false;
}

bool AdsroptionSimpleCubic4sMulti::mf_westUpTest(){
    //Check if they have only Cu
    if ( m_site->getNeighPosition(Site::NORTH)->getLabel() != "Cu") return false;
    if ( m_site->getNeighPosition(Site::WEST)->getLabel() != "Cu") return false;
    if ( m_site->getNeighPosition(Site::NORTH)->getNeighPosition(Site::WEST)->getLabel() != "Cu") return false;

 /*   cout << "T " << m_site->getID() << " ";
    cout << "N " << m_site->getNeighPosition(Site::NORTH)->getID() << " ";
    cout << "W " << m_site->getNeighPosition(Site::WEST)->getID() << " ";
    cout << "S " << m_site->getNeighPosition(Site::SOUTH)->getID() << " ";
    cout << "E " << m_site->getNeighPosition(Site::EAST)->getID() << " ";
    cout << "NW " << m_site->getNeighPosition(Site::NORTH)->getNeighPosition(Site::WEST)->getID() << " ";
    cout << "NE " << m_site->getNeighPosition(Site::NORTH)->getNeighPosition(Site::EAST)->getID() << " ";
    cout << "SW " << m_site->getNeighPosition(Site::SOUTH)->getNeighPosition(Site::WEST)->getID() << " ";
    cout << "SE " << m_site->getNeighPosition(Site::SOUTH)->getNeighPosition(Site::EAST)->getID() << " ";
    cout << endl;*/

    if ( m_site->getHeight() != m_site->getNeighPosition(Site::NORTH)->getHeight() ) return false;
    if ( m_site->getHeight() != m_site->getNeighPosition(Site::WEST)->getHeight()) return false;
    if ( m_site->getHeight() != m_site->getNeighPosition(Site::NORTH)->getNeighPosition(Site::WEST)->getHeight() ) return false;

    return true;
}

bool AdsroptionSimpleCubic4sMulti::mf_westDownTest(){
    //Check if they have only Cu
    if ( m_site->getNeighPosition(Site::SOUTH)->getLabel() != "Cu") return false;
    if ( m_site->getNeighPosition(Site::WEST)->getLabel() != "Cu") return false;
    if ( m_site->getNeighPosition(Site::SOUTH)->getNeighPosition(Site::WEST)->getLabel() != "Cu") return false;

    if ( m_site->getHeight() != m_site->getNeighPosition(Site::SOUTH)->getHeight() ) return false;
    if ( m_site->getHeight() != m_site->getNeighPosition(Site::WEST)->getHeight()) return false;
    if ( m_site->getHeight() != m_site->getNeighPosition(Site::SOUTH)->getNeighPosition(Site::WEST)->getHeight() ) return false;

    return true;
}

bool AdsroptionSimpleCubic4sMulti::mf_eastUpTest(){
    //Check if they have only Cu
    if ( m_site->getNeighPosition(Site::NORTH)->getLabel() != "Cu") return false;
    if ( m_site->getNeighPosition(Site::EAST)->getLabel() != "Cu") return false;
    if ( m_site->getNeighPosition(Site::NORTH)->getNeighPosition(Site::EAST)->getLabel() != "Cu") return false;

    if ( m_site->getHeight() != m_site->getNeighPosition(Site::NORTH)->getHeight() ) return false;
    if ( m_site->getHeight() != m_site->getNeighPosition(Site::EAST)->getHeight()) return false;
    if ( m_site->getHeight() != m_site->getNeighPosition(Site::NORTH)->getNeighPosition(Site::EAST)->getHeight() ) return false;

    return true;
}

bool AdsroptionSimpleCubic4sMulti::mf_eastDownTest(){
    //Check if they have only Cu
    if ( m_site->getNeighPosition(Site::SOUTH)->getLabel() != "Cu") return false;
    if ( m_site->getNeighPosition(Site::EAST)->getLabel() != "Cu") return false;
    if ( m_site->getNeighPosition(Site::SOUTH)->getNeighPosition(Site::EAST)->getLabel() != "Cu") return false;

    if ( m_site->getHeight() != m_site->getNeighPosition(Site::SOUTH)->getHeight() ) return false;
    if ( m_site->getHeight() != m_site->getNeighPosition(Site::EAST)->getHeight()) return false;
    if ( m_site->getHeight() != m_site->getNeighPosition(Site::SOUTH)->getNeighPosition(Site::EAST)->getHeight() ) return false;

    return true;
}

void AdsroptionSimpleCubic4sMulti::perform( Site* s )
{
    m_seAffectedSites.clear();
    m_site = s;

    //How many 4sites are availble?    //Check all and pick randomly one
    vector<string> toPerform;
    if ( mf_westUpTest() )
        toPerform.push_back( "wu");

    if ( mf_westDownTest() )
        toPerform.push_back( "wd");

    if ( mf_eastUpTest() )
        toPerform.push_back( "eu");

    if ( mf_eastDownTest() )
        toPerform.push_back( "ed");

    string p;
    if ( toPerform.size() > 1 )
        p = toPerform[ rand() % toPerform.size() ];
    else if ( toPerform.size() == 1 )
        p = toPerform[ 0 ];
    else if ( toPerform.empty() ) {
        cout << "ERROR something went wrong on adsoprtion " << endl;
        exit(0);
    }

    if ( p == "wu" )
         mf_performWU();
    else if ( p == "wd" )
        mf_performWD();
    else if ( p == "eu" )
        mf_perforEU();
    else if ( p == "ed" )
        mf_performED();
}

void AdsroptionSimpleCubic4sMulti::mf_performWU()
{
    cout << "WU" << endl;

    m_site->setLabel("CuAMD");
    m_site->getNeighPosition(Site::NORTH)->setLabel("CuAMD");
    m_site->getNeighPosition(Site::NORTH)->getNeighPosition(Site::WEST)->setLabel("AMD");
    m_site->getNeighPosition(Site::WEST)->setLabel("AMD");

    m_site->setCoupledSite( m_site->getNeighPosition(Site::NORTH) );
    m_site->getNeighPosition(Site::NORTH)->setCoupledSite( m_site );
    m_site->getNeighPosition(Site::WEST)->setCoupledSite( m_site->getNeighPosition(Site::NORTH)->getNeighPosition(Site::WEST) );
    m_site->getNeighPosition(Site::NORTH)->getNeighPosition(Site::WEST)->setCoupledSite( m_site->getNeighPosition(Site::WEST) );

    m_seAffectedSites.insert(  m_site );
    m_seAffectedSites.insert(  m_site->getNeighPosition(Site::NORTH) );
    m_seAffectedSites.insert(  m_site->getNeighPosition(Site::WEST) );
    m_seAffectedSites.insert(  m_site->getNeighPosition(Site::EAST) );
    m_seAffectedSites.insert(  m_site->getNeighPosition(Site::SOUTH) );

    m_seAffectedSites.insert(  m_site->getNeighPosition(Site::NORTH)->getNeighPosition(Site::WEST) );
    m_seAffectedSites.insert(  m_site->getNeighPosition(Site::NORTH)->getNeighPosition(Site::NORTH) );
    m_seAffectedSites.insert(  m_site->getNeighPosition(Site::NORTH)->getNeighPosition(Site::EAST) );

    m_seAffectedSites.insert(  m_site->getNeighPosition(Site::NORTH)->getNeighPosition(Site::WEST)->getNeighPosition(Site::WEST) );
    m_seAffectedSites.insert(  m_site->getNeighPosition(Site::NORTH)->getNeighPosition(Site::NORTH)->getNeighPosition(Site::WEST) );
    m_seAffectedSites.insert(  m_site->getNeighPosition(Site::NORTH)->getNeighPosition(Site::NORTH)->getNeighPosition(Site::EAST) );
    m_seAffectedSites.insert(  m_site->getNeighPosition(Site::NORTH)->getNeighPosition(Site::NORTH)->getNeighPosition(Site::WEST)->getNeighPosition(Site::WEST) );

    m_seAffectedSites.insert(  m_site->getNeighPosition(Site::WEST)->getNeighPosition(Site::SOUTH) );
    m_seAffectedSites.insert(  m_site->getNeighPosition(Site::WEST)->getNeighPosition(Site::WEST) );
    m_seAffectedSites.insert(  m_site->getNeighPosition(Site::EAST)->getNeighPosition(Site::SOUTH) );
    m_seAffectedSites.insert(  m_site->getNeighPosition(Site::WEST)->getNeighPosition(Site::WEST)->getNeighPosition(Site::SOUTH) );
}

void AdsroptionSimpleCubic4sMulti::mf_performWD()
{
    cout << "WD" << endl;

    m_site->setLabel("CuAMD");
    m_site->getNeighPosition(Site::SOUTH)->setLabel("CuAMD");
    m_site->getNeighPosition(Site::SOUTH)->getNeighPosition(Site::WEST)->setLabel("AMD");
    m_site->getNeighPosition(Site::WEST)->setLabel("AMD");

    m_site->setCoupledSite( m_site->getNeighPosition(Site::SOUTH) );
    m_site->getNeighPosition(Site::SOUTH)->setCoupledSite( m_site );
    m_site->getNeighPosition(Site::WEST)->setCoupledSite( m_site->getNeighPosition(Site::SOUTH)->getNeighPosition(Site::WEST) );
    m_site->getNeighPosition(Site::SOUTH)->getNeighPosition(Site::WEST)->setCoupledSite( m_site->getNeighPosition(Site::WEST) );

    m_seAffectedSites.insert(  m_site );
    m_seAffectedSites.insert(  m_site->getNeighPosition(Site::NORTH) );
    m_seAffectedSites.insert(  m_site->getNeighPosition(Site::WEST) );
    m_seAffectedSites.insert(  m_site->getNeighPosition(Site::EAST) );
    m_seAffectedSites.insert(  m_site->getNeighPosition(Site::SOUTH) );

    m_seAffectedSites.insert(  m_site->getNeighPosition(Site::SOUTH)->getNeighPosition(Site::WEST) );
    m_seAffectedSites.insert(  m_site->getNeighPosition(Site::NORTH)->getNeighPosition(Site::EAST) );
    m_seAffectedSites.insert(  m_site->getNeighPosition(Site::SOUTH)->getNeighPosition(Site::EAST) );
    m_seAffectedSites.insert(  m_site->getNeighPosition(Site::SOUTH)->getNeighPosition(Site::SOUTH) );
    m_seAffectedSites.insert(  m_site->getNeighPosition(Site::SOUTH)->getNeighPosition(Site::WEST)->getNeighPosition(Site::SOUTH) );
    m_seAffectedSites.insert(  m_site->getNeighPosition(Site::WEST)->getNeighPosition(Site::WEST) );
    m_seAffectedSites.insert(  m_site->getNeighPosition(Site::NORTH)->getNeighPosition(Site::WEST) );
    m_seAffectedSites.insert(  m_site->getNeighPosition(Site::NORTH)->getNeighPosition(Site::WEST)->getNeighPosition(Site::WEST) );
    m_seAffectedSites.insert(  m_site->getNeighPosition(Site::SOUTH)->getNeighPosition(Site::WEST)->getNeighPosition(Site::WEST) );
    m_seAffectedSites.insert(  m_site->getNeighPosition(Site::SOUTH)->getNeighPosition(Site::WEST)->getNeighPosition(Site::WEST)->getNeighPosition(Site::SOUTH) );
    m_seAffectedSites.insert(  m_site->getNeighPosition(Site::SOUTH)->getNeighPosition(Site::EAST)->getNeighPosition(Site::SOUTH) );
}

void AdsroptionSimpleCubic4sMulti::mf_perforEU()
{
    cout << "EU" << m_site->getID() << endl;

    m_site->setLabel("CuAMD");
    m_site->getNeighPosition(Site::NORTH)->setLabel("CuAMD");
    m_site->getNeighPosition(Site::NORTH)->getNeighPosition(Site::EAST)->setLabel("AMD");
    m_site->getNeighPosition(Site::EAST)->setLabel("AMD");

    m_site->setCoupledSite( m_site->getNeighPosition(Site::NORTH) );
    m_site->getNeighPosition(Site::NORTH)->setCoupledSite( m_site );
    m_site->getNeighPosition(Site::EAST)->setCoupledSite( m_site->getNeighPosition(Site::NORTH)->getNeighPosition(Site::EAST) );
    m_site->getNeighPosition(Site::NORTH)->getNeighPosition(Site::EAST)->setCoupledSite( m_site->getNeighPosition(Site::EAST) );


        m_seAffectedSites.insert(  m_site );
        m_seAffectedSites.insert(  m_site->getNeighPosition(Site::NORTH) );
        m_seAffectedSites.insert(  m_site->getNeighPosition(Site::WEST) );
        m_seAffectedSites.insert(  m_site->getNeighPosition(Site::EAST) );
        m_seAffectedSites.insert(  m_site->getNeighPosition(Site::SOUTH) );

        m_seAffectedSites.insert(  m_site->getNeighPosition(Site::NORTH)->getNeighPosition(Site::EAST) );
        m_seAffectedSites.insert(  m_site->getNeighPosition(Site::NORTH)->getNeighPosition(Site::EAST)->getNeighPosition(Site::EAST) );

        m_seAffectedSites.insert(  m_site->getNeighPosition(Site::EAST)->getNeighPosition(Site::EAST) );
        m_seAffectedSites.insert(  m_site->getNeighPosition(Site::EAST)->getNeighPosition(Site::EAST)->getNeighPosition(Site::SOUTH) );
        m_seAffectedSites.insert(  m_site->getNeighPosition(Site::EAST)->getNeighPosition(Site::SOUTH) );

        m_seAffectedSites.insert(  m_site->getNeighPosition(Site::NORTH)->getNeighPosition(Site::EAST)->getNeighPosition(Site::EAST)->getNeighPosition(Site::NORTH) );
        m_seAffectedSites.insert(  m_site->getNeighPosition(Site::NORTH)->getNeighPosition(Site::EAST)->getNeighPosition(Site::NORTH) );

        m_seAffectedSites.insert(  m_site->getNeighPosition(Site::NORTH)->getNeighPosition(Site::NORTH) );
        m_seAffectedSites.insert(  m_site->getNeighPosition(Site::NORTH)->getNeighPosition(Site::NORTH)->getNeighPosition(Site::WEST) );
        m_seAffectedSites.insert(  m_site->getNeighPosition(Site::NORTH)->getNeighPosition(Site::WEST) );
        m_seAffectedSites.insert(  m_site->getNeighPosition(Site::SOUTH)->getNeighPosition(Site::WEST) );


/*    m_seAffectedSites.insert(  m_site );
    m_seAffectedSites.insert(  m_site->getNeighPosition(Site::NORTH) );
    m_seAffectedSites.insert(  m_site->getNeighPosition(Site::WEST) );
    m_seAffectedSites.insert(  m_site->getNeighPosition(Site::EAST) );
    m_seAffectedSites.insert(  m_site->getNeighPosition(Site::SOUTH) );

    m_seAffectedSites.insert(  m_site->getNeighPosition(Site::NORTH)->getNeighPosition(Site::EAST) );
    m_seAffectedSites.insert(  m_site->getNeighPosition(Site::NORTH)->getNeighPosition(Site::EAST)->getNeighPosition(Site::EAST) );

    m_seAffectedSites.insert(  m_site->getNeighPosition(Site::EAST)->getNeighPosition(Site::EAST) );
    m_seAffectedSites.insert(  m_site->getNeighPosition(Site::EAST)->getNeighPosition(Site::SOUTH) );
    m_seAffectedSites.insert(  m_site->getNeighPosition(Site::NORTH)->getNeighPosition(Site::EAST)->getNeighPosition(Site::SOUTH) );

    m_seAffectedSites.insert(  m_site->getNeighPosition(Site::NORTH)->getNeighPosition(Site::EAST)->getNeighPosition(Site::EAST)->getNeighPosition(Site::NORTH) );
    m_seAffectedSites.insert(  m_site->getNeighPosition(Site::NORTH)->getNeighPosition(Site::EAST)->getNeighPosition(Site::NORTH) );

    m_seAffectedSites.insert(  m_site->getNeighPosition(Site::NORTH)->getNeighPosition(Site::NORTH) );
    m_seAffectedSites.insert(  m_site->getNeighPosition(Site::NORTH)->getNeighPosition(Site::NORTH)->getNeighPosition(Site::WEST) );
    m_seAffectedSites.insert(  m_site->getNeighPosition(Site::NORTH)->getNeighPosition(Site::WEST) );
    m_seAffectedSites.insert(  m_site->getNeighPosition(Site::SOUTH)->getNeighPosition(Site::WEST) );*/

/*    cout << "T " << m_site->getID() << endl;
    cout << "N " << m_site->getNeighPosition(Site::NORTH)->getID() << endl;
    cout << "W " << m_site->getNeighPosition(Site::WEST)->getID() << endl;
    cout << "E " << m_site->getNeighPosition(Site::EAST)->getID() << endl;
    cout << "S " << m_site->getNeighPosition(Site::SOUTH)->getID() << endl;

    cout << "NE " << m_site->getNeighPosition(Site::NORTH)->getNeighPosition(Site::EAST)->getID() << endl;
    cout << "NEE " <<m_site->getNeighPosition(Site::NORTH)->getNeighPosition(Site::EAST)->getNeighPosition(Site::EAST)->getID() << endl;

    cout << "EE " << m_site->getNeighPosition(Site::EAST)->getNeighPosition(Site::EAST)->getID() << endl;
    cout << "ES " << m_site->getNeighPosition(Site::EAST)->getNeighPosition(Site::SOUTH)->getID() << endl;
    cout << "NES " << m_site->getNeighPosition(Site::NORTH)->getNeighPosition(Site::EAST)->getNeighPosition(Site::SOUTH)->getID()  << endl;

    cout << "NEEN " << m_site->getNeighPosition(Site::NORTH)->getNeighPosition(Site::EAST)->getNeighPosition(Site::EAST)->getNeighPosition(Site::NORTH)->getID() << endl;
    cout << "NEN " << m_site->getNeighPosition(Site::NORTH)->getNeighPosition(Site::EAST)->getNeighPosition(Site::NORTH)->getID() << endl;

    cout << "NN " << m_site->getNeighPosition(Site::NORTH)->getNeighPosition(Site::NORTH)->getID() << endl;
    cout << "NNW " <<  m_site->getNeighPosition(Site::NORTH)->getNeighPosition(Site::NORTH)->getNeighPosition(Site::WEST)->getID() << endl;
    cout << "NW " << m_site->getNeighPosition(Site::NORTH)->getNeighPosition(Site::WEST)->getID() << endl;
    cout << "SW " << m_site->getNeighPosition(Site::SOUTH)->getNeighPosition(Site::WEST)->getID()  << endl;*/

}

void AdsroptionSimpleCubic4sMulti::mf_performED()
{
    cout << "ED " << m_site->getID() << endl;
    m_site->setLabel("CuAMD");
    m_site->getNeighPosition(Site::SOUTH)->setLabel("CuAMD");
    m_site->getNeighPosition(Site::SOUTH)->getNeighPosition(Site::EAST)->setLabel("AMD");
    m_site->getNeighPosition(Site::EAST)->setLabel("AMD");

    m_site->setCoupledSite( m_site->getNeighPosition(Site::SOUTH) );
    m_site->getNeighPosition(Site::SOUTH)->setCoupledSite( m_site );
    m_site->getNeighPosition(Site::EAST)->setCoupledSite( m_site->getNeighPosition(Site::SOUTH)->getNeighPosition(Site::EAST) );
    m_site->getNeighPosition(Site::SOUTH)->getNeighPosition(Site::EAST)->setCoupledSite( m_site->getNeighPosition(Site::EAST) );

    m_seAffectedSites.insert(  m_site );
    m_seAffectedSites.insert(  m_site->getNeighPosition(Site::NORTH) );
    m_seAffectedSites.insert(  m_site->getNeighPosition(Site::WEST) );
    m_seAffectedSites.insert(  m_site->getNeighPosition(Site::EAST) );
    m_seAffectedSites.insert(  m_site->getNeighPosition(Site::SOUTH) );

    m_seAffectedSites.insert(  m_site->getNeighPosition(Site::NORTH)->getNeighPosition(Site::WEST) );
    m_seAffectedSites.insert(  m_site->getNeighPosition(Site::NORTH)->getNeighPosition(Site::EAST) );
    m_seAffectedSites.insert(  m_site->getNeighPosition(Site::NORTH)->getNeighPosition(Site::EAST)->getNeighPosition(Site::EAST) );

    m_seAffectedSites.insert(  m_site->getNeighPosition(Site::EAST)->getNeighPosition(Site::SOUTH) );
    m_seAffectedSites.insert(  m_site->getNeighPosition(Site::EAST)->getNeighPosition(Site::EAST) );
    m_seAffectedSites.insert(  m_site->getNeighPosition(Site::EAST)->getNeighPosition(Site::EAST)->getNeighPosition(Site::SOUTH) );
    m_seAffectedSites.insert(  m_site->getNeighPosition(Site::EAST)->getNeighPosition(Site::EAST)->getNeighPosition(Site::SOUTH)->getNeighPosition(Site::SOUTH) );

    m_seAffectedSites.insert(  m_site->getNeighPosition(Site::EAST)->getNeighPosition(Site::SOUTH)->getNeighPosition(Site::SOUTH) );
    m_seAffectedSites.insert(  m_site->getNeighPosition(Site::SOUTH)->getNeighPosition(Site::SOUTH) );
    m_seAffectedSites.insert(  m_site->getNeighPosition(Site::SOUTH)->getNeighPosition(Site::SOUTH)->getNeighPosition(Site::WEST) );
    m_seAffectedSites.insert(  m_site->getNeighPosition(Site::SOUTH)->getNeighPosition(Site::WEST) );


/*    cout << "T " << m_site->getID() << endl;
    cout << "N " << m_site->getNeighPosition(Site::NORTH)->getID() << endl;

    cout << "W " << m_site->getNeighPosition(Site::WEST)->getID() << endl;
    cout << "E " << m_site->getNeighPosition(Site::EAST)->getID() << endl;
    cout << "S " << m_site->getNeighPosition(Site::SOUTH)->getID() << endl;

    cout << "NW " << m_site->getNeighPosition(Site::NORTH)->getNeighPosition(Site::WEST)->getID() << endl;
    cout << "NE " << m_site->getNeighPosition(Site::NORTH)->getNeighPosition(Site::EAST)->getID() << endl;
    cout << "NEE " << m_site->getNeighPosition(Site::NORTH)->getNeighPosition(Site::EAST)->getNeighPosition(Site::EAST) ->getID() << endl;

    cout << "ES " << m_site->getNeighPosition(Site::EAST)->getNeighPosition(Site::SOUTH)->getID() << endl;
    cout << "EE " << m_site->getNeighPosition(Site::EAST)->getNeighPosition(Site::EAST)->getID()  << endl;
    cout << "EES " << m_site->getNeighPosition(Site::EAST)->getNeighPosition(Site::EAST)->getNeighPosition(Site::SOUTH)->getID() << endl;
    cout << "EESS " << m_site->getNeighPosition(Site::EAST)->getNeighPosition(Site::EAST)->getNeighPosition(Site::SOUTH)->getNeighPosition(Site::SOUTH)->getID() << endl;

    cout << "ESS " << m_site->getNeighPosition(Site::EAST)->getNeighPosition(Site::SOUTH)->getNeighPosition(Site::SOUTH)->getID() << endl;
    cout << "SS " << m_site->getNeighPosition(Site::SOUTH)->getNeighPosition(Site::SOUTH)->getID() << endl;
    cout << "SSW " <<m_site->getNeighPosition(Site::SOUTH)->getNeighPosition(Site::SOUTH)->getNeighPosition(Site::WEST)->getID() << endl;
    cout << "SW " << m_site->getNeighPosition(Site::SOUTH)->getNeighPosition(Site::WEST)->getID()  << endl;*/


}

double AdsroptionSimpleCubic4sMulti::getProbability(){

    //These must trenafered in the global definitions
    double Na = 6.0221417930e+23;		// Avogadro's number [1/mol]
    double P = any_cast<double>(m_mParams["P"]);					// [Pa]
    double T = any_cast<double>(m_mParams["T"]); //500;						// [K]
    double k = any_cast<double>(m_mParams["k"]); // 1.3806503e-23;			// Boltzmann's constant [j/K]
    double s0 = any_cast<double>(m_mParams["s0"]); //0.1;
    double C_tot = any_cast<double>(m_mParams["C_tot"]);			// [sites/m^2] Vlachos code says [moles sites/m^2]
    // Ctot for copper: 2e+19
    double m = (409.1e-3)/Na;				// [kg/mol] this is the molecular weight
    double y = any_cast<double>(m_mParams["f"]);					// Mole fraction of the precursor on the wafer

    return s0*y*P/(C_tot*sqrt(2.0e0*3.14159265*m*k*T) );
}

}
