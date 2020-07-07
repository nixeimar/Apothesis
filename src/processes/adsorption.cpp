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

#include "adsorption.h"
#include "register.cpp"
#include "parameters.h"

namespace MicroProcesses{

Adsorption::Adsorption
(
  string species,
  double stickingCoeffs,
  double massFraction
)
:
m_sName("Adsorption"),
m_iNeighNum(0), 
m_apothesis(0),
m_adsorptionSpecies(species),
m_stickingCoeffs(stickingCoeffs),
m_massfraction(massFraction),
m_canDesorb(false)
{
;  
}

Adsorption::~Adsorption(){}

string Adsorption::getName(){ return m_sName; }

string Adsorption::getSpecies(){ return m_adsorptionSpecies; }

//This should be called only once in the initialization
void Adsorption::activeSites(Lattice* lattice)
  {
  m_pLattice = lattice;
  vector< Site* > vSites = m_pLattice->getSites();

  for ( int i = 0; i < m_pLattice->getSize(); i++)
    if ( vSites[ i ]->getID()%2 != 0 ) {
      m_lAdsSites.push_back( vSites[ i ] );
      vSites[ i ]->addProcess( this );
      }
  }

void Adsorption::selectSite()
{
  /* This comes from random i.e. picking from the available list for adsorption randomly */
  int y = rand()%getActiveList().size();
  int counter = 0;
  list<Site* >::iterator site = m_lAdsSites.begin();
  for (; site != m_lAdsSites.end(); site++ ){
    if ( counter == y )
      m_site = (*site);
    counter++;
    }
}

void Adsorption::setProcessMap( map< Process*, list<Site* >* >* ){}

void Adsorption::perform()
{
  selectSite();

  int height = m_site->getHeight();
  height = height + 2;
  m_site->setHeight( height);
  // How to set a random adsorption species? Get a pointer to apothesis in order to find/get access to the species list
  //m_site->setSpecies(m_adsorptionSpecies[0]);
  
  mf_removeFromList();
  mf_updateNeighNum();

  // Add desorption site to Desorption class
  if (canDesorb())
  {
    getDesorption()->mf_addToList(m_site);
  }

  /// Check if there are available sites that it can be performed
  if (m_lAdsSites.size() == 0)
  {
    cout << "No more "<<getName()<< " site is available. Exiting..." << endl;
    m_apothesis->pErrorHandler->error_simple_msg( "No "+ getName() + " site is available.");
    EXIT;
  }
}

void Adsorption::mf_removeFromList() { m_lAdsSites.remove( m_site); m_site->removeProcess( this ); }

void Adsorption::mf_addToList(Site *s) { m_lAdsSites.push_back( s); }

void Adsorption::mf_updateNeighNum()
{
  int siteHeight = m_site->getHeight();

  bool isActiveEAST = false;
  isActiveEAST = ( siteHeight == m_site->getNeighPosition( Site::EAST )->getHeight()  && \
                   siteHeight == m_site->getNeighPosition( Site::EAST_DOWN)->getHeight() && \
                   siteHeight == m_site->getNeighPosition( Site::EAST_UP)->getHeight() );

  bool isActiveWEST = false;
  isActiveWEST = ( siteHeight == m_site->getNeighPosition( Site::WEST )->getHeight() && \
                   siteHeight == m_site->getNeighPosition( Site::WEST_DOWN)->getHeight() && \
                   siteHeight == m_site->getNeighPosition( Site::WEST_UP)->getHeight() );

  bool isActiveNORTH = false;
  isActiveNORTH = ( siteHeight == m_site->getNeighPosition( Site::WEST_UP )->getHeight() && \
                    siteHeight == m_site->getNeighPosition( Site::EAST_UP)->getHeight() && \
                    siteHeight == m_site->getNeighPosition( Site::NORTH)->getHeight() );

  bool isActiveSOUTH = false;
  isActiveSOUTH = ( siteHeight == m_site->getNeighPosition( Site::WEST_DOWN )->getHeight() && \
                    siteHeight == m_site->getNeighPosition( Site::EAST_DOWN)->getHeight() &&  \
                    siteHeight == m_site->getNeighPosition( Site::SOUTH)->getHeight() );

  // Store the activated sites
  if ( isActiveEAST )
    mf_addToList( m_site->getActivationSite( Site::ACTV_EAST ));

  if ( isActiveWEST )
    mf_addToList( m_site->getActivationSite( Site::ACTV_WEST ));

  if ( isActiveNORTH )
    mf_addToList( m_site->getActivationSite( Site::ACTV_NORTH ));

  if ( isActiveSOUTH )
    mf_addToList( m_site->getActivationSite( Site::ACTV_SOUTH ));
}

const double Adsorption::getMassFraction()
{
  return m_massfraction;
}

list<Site*> Adsorption::getActiveList()
{
  return m_lAdsSites;
}

void Adsorption::test()
{
  cout << m_lAdsSites.size() << endl;
}

double Adsorption::getProbability()
{
  /* These are parameters values (I/O) */
  double dNavogadro = m_apothesis->pParameters->dAvogadroNum;
  double dPres = m_apothesis->pParameters->getPressure();
  double dTemp = m_apothesis->pParameters->getTemperature();
  double dkBoltz = m_apothesis->pParameters->dkBoltz;
  
  int index = m_apothesis->findSpeciesIndex(m_adsorptionSpecies);

  double dmass = 27e-3/dNavogadro;
  double dpi = 3.14159265;
  double dstick = m_stickingCoeffs;
  double dCites = 1.4e+19;
  double dy = getMassFraction();

  /* Adsorption probability see Lam and Vlachos */
  double dflux = dstick*dPres*dy/(dCites*sqrt(2.0*dpi*dmass*dkBoltz*dTemp));

  if ( m_lAdsSites.size() !=0 )
    return m_lAdsSites.size()*dflux;
  else
    return 0.0;
}

Desorption* Adsorption::getDesorption()
{
  return m_pDesorption;
}

void Adsorption::setDesorptionPointer(Desorption* d)
{
  m_pDesorption = d;
}

bool Adsorption::canDesorb()
{
  return m_canDesorb;
}

// Set desorption boolean variable to be true 
void Adsorption::setDesorption()
{
  m_canDesorb = true; 
}

}
