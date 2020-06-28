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
m_desorptionEnergy(energy),
m_desorptionFrequency(frequency)
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
  // Find the species on the site

  /* These are parameters values (I/O) */
  double dTemp = m_apothesis->pParameters->getTemperature();
  double dkBoltz = m_apothesis->pParameters->dkBoltz;

  // how to find the number of neighbours?
 //   {
 //     cout << "No more "<<m_vProcesses[0]->getName()<< " site is available. Exiting..." << endl;
 //     pErrorHandler->error_simple_msg( "No "+ m_vProcesses[0]->getName() + " site is available. Last time step: " + to_string( i ) );
 //     EXIT;
 //   }

  double n = 1;
  int index = -1;
  
  // Find index of individual species
  for (int i = 0; i < m_desorptionSpecies.size(); i++)
  {
    //if(!m_desorptionSpecies[i].compare(species))
    //{
    //  index = i; 
    //}
  }

  double freq = m_desorptionFrequency[index];
  double energy = m_desorptionEnergy[index];

  // TODO: Error message if index = -1

  /* Desorption probability see Lam and Vlachos  */
  double dflux = freq*exp(-n*energy/(dkBoltz*dTemp));

  if ( m_lAdsSites.size() !=0 )
    return m_lAdsSites.size()*dflux;
  else
    return 0.0;
}

const vector<string> Desorption::getDesorptionSpecies()
{
  return m_desorptionSpecies;
}

const vector<double> Desorption::getDesorptionEnergy()
{
  return m_desorptionEnergy;
}

const vector<double> Desorption::getDesorptionFrequency()
{
  return m_desorptionFrequency;
}

}
