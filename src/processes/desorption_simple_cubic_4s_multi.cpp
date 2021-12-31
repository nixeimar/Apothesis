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


#include "desorption_simple_cubic_4s_multi.h"

namespace MicroProcesses
{

REGISTER_PROCESS_IMPL( DesorptionSimpleCubic4sMulti )

DesorptionSimpleCubic4sMulti::DesorptionSimpleCubic4sMulti():m_site(0)
{
    /* initialize random seed: */
   // srand(2000);
    srand (time(NULL));
}

DesorptionSimpleCubic4sMulti::~DesorptionSimpleCubic4sMulti(){}

bool DesorptionSimpleCubic4sMulti::rules( Site* s )
{
    m_site = s;
    if ( !s->getCoupledSite() ) return false;
    if ( mf_calculateNeighs() != m_iNeighs ) return false;

    return true;
}

int DesorptionSimpleCubic4sMulti::mf_calculateNeighs()
{
    int iCountNeighs = 0;

    if ( m_site->getNeighPosition(Site::EAST)->getHeight() >= m_site->getHeight() + 1 &&
         (m_site->getNeighPosition(Site::EAST)->getLabel() == "Cu" || m_site->getNeighPosition(Site::EAST)->getLabel() == "CuAMD" ) )
        iCountNeighs++;

    if ( m_site->getNeighPosition(Site::WEST)->getHeight() >= m_site->getHeight() + 1 &&
         (m_site->getNeighPosition(Site::WEST)->getLabel() == "Cu" || m_site->getNeighPosition(Site::WEST)->getLabel() == "CuAMD" ) )
        iCountNeighs++;

    if ( m_site->getCoupledSite()->getNeighPosition(Site::EAST)->getHeight() >= m_site->getCoupledSite()->getHeight() + 1 &&
         (m_site->getCoupledSite()->getNeighPosition(Site::EAST)->getLabel() == "Cu" || m_site->getCoupledSite()->getNeighPosition(Site::EAST)->getLabel() == "CuAMD" ) )
        iCountNeighs++;

    if ( m_site->getCoupledSite()->getNeighPosition(Site::WEST)->getHeight() >= m_site->getCoupledSite()->getHeight() + 1
         && (m_site->getCoupledSite()->getNeighPosition(Site::WEST)->getLabel() == "Cu" || m_site->getCoupledSite()->getNeighPosition(Site::WEST)->getLabel() == "CuAMD" ) )
        iCountNeighs++;

    if ( m_site->getCoupledSite() == m_site->getNeighPosition(Site::SOUTH) ) {
        if ( m_site->getNeighPosition(Site::NORTH)->getHeight() >= m_site->getHeight() + 1 &&
             (m_site->getNeighPosition(Site::NORTH)->getLabel() == "Cu" || m_site->getNeighPosition(Site::NORTH)->getLabel() == "CuAMD" )  )
            iCountNeighs++;

        if ( m_site->getCoupledSite()->getNeighPosition(Site::SOUTH)->getNeighPosition(Site::SOUTH)->getHeight() >= m_site->getCoupledSite()->getHeight() + 1 &&
             (m_site->getCoupledSite()->getNeighPosition(Site::SOUTH)->getLabel() == "Cu" ||
                 m_site->getCoupledSite()->getNeighPosition(Site::SOUTH)->getLabel() == "CuAMD" ) )
            iCountNeighs++;
    }
    else if ( m_site->getCoupledSite() == m_site->getNeighPosition(Site::NORTH) ) {
        if ( m_site->getNeighPosition(Site::SOUTH)->getHeight() >= m_site->getHeight() + 1
             && ( m_site->getNeighPosition(Site::SOUTH)->getLabel() == "Cu" ||m_site->getNeighPosition(Site::SOUTH)->getLabel() == "CuAMD" ))
            iCountNeighs++;

        if ( m_site->getCoupledSite()->getNeighPosition(Site::NORTH)->getNeighPosition(Site::NORTH)->getHeight() >= m_site->getCoupledSite()->getHeight() + 1 &&
             ( m_site->getCoupledSite()->getNeighPosition(Site::NORTH)->getNeighPosition(Site::NORTH)->getLabel() == "Cu" ||
                   m_site->getCoupledSite()->getNeighPosition(Site::NORTH)->getNeighPosition(Site::NORTH)->getLabel() == "CuAMD" ) )
            iCountNeighs++;
    }

    return iCountNeighs;
}

void DesorptionSimpleCubic4sMulti::perform( Site* s )
{
    m_site = s;

    if ( !s->getCoupledSite() ){
        cout << "Oups something is Fucked up ..." << endl;
        exit(0);
    }

    m_seAffectedSites.clear();
    m_seAffectedSites.insert( s );
    m_seAffectedSites.insert( s->getNeighPosition(Site::NORTH) );
    m_seAffectedSites.insert( s->getNeighPosition(Site::EAST) );
    m_seAffectedSites.insert( s->getNeighPosition(Site::WEST) );
    m_seAffectedSites.insert( s->getNeighPosition(Site::SOUTH) );
    m_seAffectedSites.insert( s->getNeighPosition(Site::NORTH)->getNeighPosition(Site::EAST) );
    m_seAffectedSites.insert( s->getNeighPosition(Site::NORTH)->getNeighPosition(Site::WEST) );
    m_seAffectedSites.insert( s->getNeighPosition(Site::SOUTH)->getNeighPosition(Site::EAST) );
    m_seAffectedSites.insert( s->getNeighPosition(Site::SOUTH)->getNeighPosition(Site::WEST) );

    m_seAffectedSites.insert( s->getCoupledSite() );
    m_seAffectedSites.insert( s->getCoupledSite()->getNeighPosition(Site::NORTH) );
    m_seAffectedSites.insert( s->getCoupledSite()->getNeighPosition(Site::EAST) );
    m_seAffectedSites.insert( s->getCoupledSite()->getNeighPosition(Site::WEST) );
    m_seAffectedSites.insert( s->getCoupledSite()->getNeighPosition(Site::SOUTH) );
    m_seAffectedSites.insert( s->getCoupledSite()->getNeighPosition(Site::NORTH)->getNeighPosition(Site::EAST) );
    m_seAffectedSites.insert( s->getCoupledSite()->getNeighPosition(Site::NORTH)->getNeighPosition(Site::WEST) );
    m_seAffectedSites.insert( s->getCoupledSite()->getNeighPosition(Site::SOUTH)->getNeighPosition(Site::EAST) );
    m_seAffectedSites.insert( s->getCoupledSite()->getNeighPosition(Site::SOUTH)->getNeighPosition(Site::WEST) );

    if (s->getLabel() == "CuAMD"){
        s->getCoupledSite()->increaseHeight(1);
        s->getCoupledSite()->setLabel("Cu");

        s->increaseHeight(1);
        s->setLabel("Cu");
    }
    else if (s->getLabel() == "AMD") {
        s->getCoupledSite()->setLabel("Cu");
        s->setLabel("Cu");
    }

    s->getCoupledSite()->setCoupledSite( 0 );
    s->setCoupledSite( 0 );
}

double DesorptionSimpleCubic4sMulti::getProbability(){

    double Na = 6.0221417930e+23;		// Avogadro's number [1/mol]
    //double P = any_cast<double>(m_mParams["P"]);					// [Pa]
    double T = any_cast<double>(m_mParams["T"]); //500;						// [K]
    double k = any_cast<double>(m_mParams["k"]); // 1.3806503e-23;			// Boltzmann's constant [j/K]
    double s0 = any_cast<double>(m_mParams["s0"]); //0.1;
    double C_tot = any_cast<double>(m_mParams["C_tot"]);			// [sites/m^2] Vlachos code says [moles sites/m^2]
    // Ctot for copper: 2e+19
    double m = (141.094e-3)/Na;				// [kg/mol] this is the molecular weight
    double y = 0.0000593894333333333; //any_cast<double>(m_mParams["f"]);					// Mole fraction of the precursor on the wafer
    double h = 6.62607004e-34; //m2 kg / s
    cout <<  (k*T/h)*exp(-(m_iNeighs*40000.)/(Na*k*T)) << endl; // this is 1.05968

    return  (k*T/h)*exp(-(m_iNeighs*40000.)/(Na*k*T));
    //For IKY: return 0.03;
//    return 150000./40.;// 1.5; //3.0; //s0*y*P/(C_tot*sqrt(2.0e0*3.14159265*m*k*T) );
}

}
