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
    Apothesis* instance,
    string speciesName,
    Species* species,
    double energy,
    double frequency
)
:
m_sName("Desorption"),
m_iNeighNum(0), 
m_apothesis(instance),
m_desorptionSpeciesName(speciesName),
m_desorptionSpecies(species),
m_desorptionEnergy(energy),
m_desorptionFrequency(frequency),
m_maxNeighbours(5)
{
  // Initialize list to state number of sites with n neighbours
  for(int i = 0; i < m_maxNeighbours; ++i)
  {
    m_numNeighbours.push_back(0);
  }

  m_probabilities = generateProbabilities();
}

Desorption::~Desorption(){}

string Desorption::getName(){ return m_sName; }

const string Desorption::getSpeciesName(){ return m_desorptionSpeciesName; }

const Species* Desorption::getSpecies(){ return m_desorptionSpecies; }

//This should be called only once in the initialization
void Desorption::activeSites(Lattice* lattice)
{

  m_pLattice = lattice;
  vector< Site* > vSites = m_pLattice->getSites();

  for ( int i = 0; i < m_pLattice->getSize(); i++)
    if ( vSites[ i ]->getID()%2 != 0 ) {
      //m_lDesSites.push_back( vSites[ i ] );
      vSites[ i ]->addProcess( this );
      }
}

void Desorption::selectSite()
{
  /* This comes from random i.e. picking from the available list for Desorption randomly */
  int len = getActiveList().size();
  int y = rand()%len;
  int counter = 0;
  list<Site* >::iterator site = m_lDesSites.begin();
  for (; site != m_lDesSites.end(); site++ ){
    if ( counter == y )
      m_site = (*site);
    counter++;
    }
}

void Desorption::setProcessMap( map< Process*, list<Site* >* >* ){}

void Desorption::perform()
{
  int height = m_site->getHeight();
  height = height - 2;
  m_site->setHeight( height);

  int numNeighbours = m_site->getNeighboursNum();
  mf_removeFromList();
  mf_updateNeighNum();

  // Remove count from list of sites that have n number of neighbours
  updateSiteCounter(numNeighbours, false);
}

void Desorption::mf_removeFromList() 
{ 
  m_lDesSites.remove( m_site); 
  //TODO: Is this necessary?
  m_site->removeProcess( this ); 
}

void Desorption::mf_addToList(Site *s) { m_lDesSites.push_back( s); }

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
  return m_lDesSites;
}

void Desorption::test()
{
  cout << m_lDesSites.size() << endl;
}


double Desorption::getProbability()
{
  if ( m_lDesSites.size() == 0 )
  {
    return 0;
  }
  
  double prob = 0;

  // Calculate probability for each possible value of n
  for (int i = 0; i < m_maxNeighbours; ++i)
  {
    prob += m_probabilities[i] * m_numNeighbours[i];
  }

  return prob;
}

vector<double> Desorption::generateProbabilities()
{
  /* These are parameters values (I/O) */
  double dTemp = m_apothesis->pParameters->getTemperature();
  double dkBoltz = m_apothesis->pParameters->dkBoltz;

  double freq = m_desorptionFrequency;
  double energy = m_desorptionEnergy;

  vector<double> prob;
  /* Desorption probability see Lam and Vlachos  */
  for (int n = 1; n <= m_maxNeighbours; ++n)
  {
    prob.push_back( freq*exp(-n*energy/(dkBoltz*dTemp)));
  }

  return prob;
}

const double Desorption::getDesorptionEnergy()
{
  return m_desorptionEnergy;
}

const double Desorption::getDesorptionFrequency()
{
  return m_desorptionFrequency;
}

Adsorption* Desorption::getAdsorption()
{
  return m_pAdsorption;
}

void Desorption::setAdsorptionPointer(Adsorption* a)
{
  m_pAdsorption = a;
}

void Desorption::updateSiteCounter(int neighbours, bool addOrRemove)
{
  // Updates list of number of neighbours each possible site has
  if (addOrRemove)
  {
    m_numNeighbours[neighbours-1]++;
  }
  else
  {
    if (m_numNeighbours[neighbours-1] > 0)
    {
      m_numNeighbours[neighbours-1]--;
    }
  }
}

void Desorption::setSite(Site* s)
{
  m_site = s;
}

}
