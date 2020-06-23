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

#include "desorption.h"
#include "register.cpp"
#include "parameters.h"

namespace MicroProcesses{

Desorption::Desorption
(
    vector<string> species,
    vector<double> energy,
    vector<double> frequency
)
:
m_sName("Desorption"),
m_iNeighNum(0), 
m_apothesis(0),
m_desorptionSpecies(species),
m_energy(energy),
m_frequency(frequency)
{
;  
}

Desorption::~Desorption(){}

string Desorption::getName(){ return m_sName; }

//This should be called only once in the initialization
void Desorption::activeSites(Lattice* lattice)
  {
  m_pLattice = lattice;
  vector< Site* > vSites = m_pLattice->getSites();

  for ( int i = 0; i < m_pLattice->getSize(); i++)
    if ( vSites[ i ]->getID()%2 != 0 ) {
      m_lAdsSites.push_back( vSites[ i ] );
      vSites[ i ]->addProcess( this );
      }
  }

void Desorption::selectSite()
{
  /* This comes from random i.e. picking from the available list for Desorption randomly */
  int y = rand()%getActiveList().size();
  int counter = 0;
  list<Site* >::iterator site = m_lAdsSites.begin();
  for (; site != m_lAdsSites.end(); site++ ){
    if ( counter == y )
      m_site = (*site);
    counter++;
    }
}

void Desorption::setProcessMap( map< Process*, list<Site* >* >* ){}

void Desorption::perform()
{
  int height = m_site->getHeight();
  height = height + 2;
  m_site->setHeight( height);
  mf_removeFromList();
  mf_updateNeighNum();
}

void Desorption::mf_removeFromList() { m_lAdsSites.remove( m_site); m_site->removeProcess( this ); }

void Desorption::mf_addToList(Site *s) { m_lAdsSites.push_back( s); }

void Desorption::mf_updateNeighNum()
{
  bool isActiveEAST = false;
  isActiveEAST = ( m_site->getHeight() == m_site->getNeighPosition( Site::EAST )->getHeight()  && \
                   m_site->getHeight() == m_site->getNeighPosition( Site::EAST_DOWN)->getHeight() && \
                   m_site->getHeight() == m_site->getNeighPosition( Site::EAST_UP)->getHeight() );

  bool isActiveWEST = false;
  isActiveWEST = ( m_site->getHeight() == m_site->getNeighPosition( Site::WEST )->getHeight() && \
                   m_site->getHeight() == m_site->getNeighPosition( Site::WEST_DOWN)->getHeight() && \
                   m_site->getHeight() == m_site->getNeighPosition( Site::WEST_UP)->getHeight() );

  bool isActiveNORTH = false;
  isActiveNORTH = ( m_site->getHeight() == m_site->getNeighPosition( Site::WEST_UP )->getHeight() && \
                    m_site->getHeight() == m_site->getNeighPosition( Site::EAST_UP)->getHeight() && \
                    m_site->getHeight() == m_site->getNeighPosition( Site::NORTH)->getHeight() );

  bool isActiveSOUTH = false;
  isActiveSOUTH = ( m_site->getHeight() == m_site->getNeighPosition( Site::WEST_DOWN )->getHeight() && \
                    m_site->getHeight() == m_site->getNeighPosition( Site::EAST_DOWN)->getHeight() &&  \
                    m_site->getHeight() == m_site->getNeighPosition( Site::SOUTH)->getHeight() );

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

list<Site*> Desorption::getActiveList()
{
  return m_lAdsSites;
}

void Desorption::test()
{
  cout << m_lAdsSites.size() << endl;
}

double Desorption::getProbability()
{
  /* These are parameters values (I/O) */
  double dNavogadro = m_apothesis->pParameters->dAvogadroNum;
  double dPres = m_apothesis->pParameters->getPressure();
  double dTemp = m_apothesis->pParameters->getTemperature();
  double dkBoltz = m_apothesis->pParameters->dkBoltz;

  double dmass = 27e-3/dNavogadro;
  double dpi = 3.14159265;
  double dstick = 0.5;
  double dCites = 1.4e+19;
  double dy = 3.05e-3;

  /* Desorption probability see Lam and Vlachos  */
  double dflux = dstick*dPres*dy/(dCites*sqrt(2.0*dpi*dmass*dkBoltz*dTemp));

  if ( m_lAdsSites.size() !=0 )
    return m_lAdsSites.size()*dflux;
  else
    return 0.0;
}

}
