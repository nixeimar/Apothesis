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

#include "diffusion.h"
#include "desorption.h"
#include "adsorption.h"
#include "register.cpp"
#include "parameters.h"
#include <cmath>
#include <algorithm>

namespace MicroProcesses{

/// Constructor
Diffusion::Diffusion
(
    Apothesis* instance,
    string species,
    double energy,
    double frequency
)
:
m_sName("Diffusion"),
m_iNeighNum(0), 
m_apothesis(instance),
m_diffusionSpecies(species),
m_diffusionEnergy(energy),
m_diffusionFrequency(frequency),
m_pDesorption(0),
m_pAdsorption(0),
m_maxNeighbours(5)
{
  for(int i = 0; i < m_maxNeighbours; ++i)
  {
    m_numNeighbours.push_back(0);
  }

  m_probabilities = generateProbabilities();
}

Diffusion::~Diffusion(){;}

void Diffusion::init(){;}

string Diffusion::getName(){ return m_sName; }

void Diffusion::activeSites( Lattice* lattice){
  m_pLattice = lattice;
  vector< Site* > vSites = m_pLattice->getSites();

  for ( int i = 0; i < m_pLattice->getSize(); i++)
    if ( vSites[ i ]->getID()%2 != 0 ) {
     // m_lDiffSites.push_back( vSites[ i ] );
      vSites[ i ]->addProcess( this );
    }
}

void Diffusion::selectSite()
{
  /* This comes from random i.e. picking from the available list for diffusion randomly */
  int y = rand()%getActiveList().size();
  int counter = 0;
  list<Site* >::iterator site = m_lDiffSites.begin();
  for (; site != m_lDiffSites.end(); site++ )
  {
    if ( counter == y )
      m_site = (*site);
    counter++;
  }
}

Site* Diffusion::chooseNeighbour(vector<Site*> neighbours)
{
  /* This comes from random i.e. picking from the available list for diffusion randomly */
  int y = rand()%neighbours.size();
  int counter = 0;
  vector<Site* >::iterator site = neighbours.begin();
  for (; site != neighbours.end(); ++site )
  {
    if ( counter == y )
      return (*site);
    counter++;
  }
}

void Diffusion::perform()
{ 
  // Choose a neighbour
  vector<Site*> neighbours = m_site->getNeighs();
  Site* diffuseTo = chooseNeighbour(neighbours);

  // Get current number of neighbours
  int currentNeighbourNum = m_site->getNeighboursNum();
  
  // Desorb, then adsorb to that neighbour
  Desorption* d = getDesorption();
  d->setSite(m_site);
  d->perform();

  Adsorption* a = getAdsorption();
  a->setSite(diffuseTo);
  a->perform();

  // update number of neighbours for all adjacent 
  vector<Site*> neigh = m_site->getNeighs();
  for(vector<Site*> :: iterator itr = neigh.begin(); itr != neigh.end(); ++itr)
  {
    Site* s = *itr;
    // TODO: how to update number of neighbours for each site
    mf_updateNeighNum(s);
  }

  // update number of possible diffusion sites
  updateSiteCounter(currentNeighbourNum, false);
  int newNeighbourNum = diffuseTo->getNeighboursNum();
  updateSiteCounter(newNeighbourNum, true);

  // get new number of neighbours
  mf_removeFromList();

  //TODO: Is this function still valid for diffusion process?
  mf_updateNeighNum(m_site);
}

void Diffusion::mf_removeFromList() { m_lDiffSites.remove( m_site); m_site->removeProcess( this ); }

void Diffusion::mf_addToList(Site *s) 
{
  //TODO: Will this be better with a hashmap/map of some kind? for now, just do a search
  if (find(m_lDiffSites.begin(), m_lDiffSites.end(),s)==m_lDiffSites.end()) 
    m_lDiffSites.push_back(s); 
}

void Diffusion::mf_updateNeighNum(Site* site)
{
  bool isActiveEAST = false;
  isActiveEAST = ( site->getHeight() == site->getNeighPosition( Site::EAST )->getHeight()  && \
                   site->getHeight() == site->getNeighPosition( Site::EAST_DOWN)->getHeight() && \
                   site->getHeight() == site->getNeighPosition( Site::EAST_UP)->getHeight() );

  bool isActiveWEST = false;
  isActiveWEST = ( site->getHeight() == site->getNeighPosition( Site::WEST )->getHeight() && \
                   site->getHeight() == site->getNeighPosition( Site::WEST_DOWN)->getHeight() && \
                   site->getHeight() == site->getNeighPosition( Site::WEST_UP)->getHeight() );

  bool isActiveNORTH = false;
  isActiveNORTH = ( site->getHeight() == site->getNeighPosition( Site::WEST_UP )->getHeight() && \
                    site->getHeight() == site->getNeighPosition( Site::EAST_UP)->getHeight() && \
                    site->getHeight() == site->getNeighPosition( Site::NORTH)->getHeight() );

  bool isActiveSOUTH = false;
  isActiveSOUTH = ( site->getHeight() == site->getNeighPosition( Site::WEST_DOWN )->getHeight() && \
                    site->getHeight() == site->getNeighPosition( Site::EAST_DOWN)->getHeight() &&  \
                    site->getHeight() == site->getNeighPosition( Site::SOUTH)->getHeight() );

  // Store the activated sites
  if ( isActiveEAST )
    mf_addToList( site->getActivationSite( Site::ACTV_EAST ));

  if ( isActiveWEST )
    mf_addToList( site->getActivationSite( Site::ACTV_WEST ));

  if ( isActiveNORTH )
    mf_addToList( site->getActivationSite( Site::ACTV_NORTH ));

  if ( isActiveSOUTH )
    mf_addToList( site->getActivationSite( Site::ACTV_SOUTH ));
}

double Diffusion::getProbability()
{
  if ( m_lDiffSites.size() == 0 )
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

list<Site* > Diffusion::getActiveList()
{
  //TODO: Figure out way to keep track of this list! from adsorption and desorption sides
  return m_lDiffSites;
}

void Diffusion::setProcessMap(map< Process*, list<Site* >* >* procMap )
{
  m_pProcessMap = procMap;
  (*m_pProcessMap)[ this] = &m_lDiffSites;
}

void Diffusion::test()
{
  cout << m_lDiffSites.size() << endl;
}

vector<double> Diffusion::generateProbabilities()
{
  /* These are parameters values (I/O) */
  double dTemp = m_apothesis->pParameters->getTemperature();
  double dkBoltz = m_apothesis->pParameters->dkBoltz;

  //these parameters are now given constant values for now.
  //later they will be taken from the input file.
  double v0 = m_diffusionFrequency;
  double E = m_diffusionEnergy;
  double Em = 0;

  vector<double> prob;
  /* Desorption probability see Lam and Vlachos  */
  for (int n = 1; n <= m_maxNeighbours; ++n)
  {
    prob.push_back(-v0*exp((E-Em)/(dkBoltz*dTemp))*exp(-n*E/(dkBoltz*dTemp)));
  }

  return prob;
}

void Diffusion::setAdsorptionPointer(Adsorption* a)
{
  m_pAdsorption = a;
}

void Diffusion::setDesorptionPointer(Desorption* d)
{
  m_pDesorption = d;
}

Adsorption* Diffusion::getAdsorption()
{
  return m_pAdsorption;
}

Desorption* Diffusion::getDesorption()
{
  return m_pDesorption;
}

void Diffusion::updateSiteCounter(int neighbours, bool addOrRemove)
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

}
