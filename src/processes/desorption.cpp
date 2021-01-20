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
m_maxNeighbours(5), //TODO: initialize maxneighbours
m_canDiffuse(false)
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
  if (m_site->getSpecies().size() == 1)
  {
    int height = m_site->getHeight();
    height = height - 2;
    m_site->setHeight( height);
  }
  
  int numNeighbours = m_site->getNeighboursNum();

  m_site->removeSpecies(m_apothesis->getSpecies(m_desorptionSpeciesName));

  
  // If there are no longer any species that can be desorbed, remove from list
  if (m_site->getSpecies().size() == 0)
  {
    mf_removeFromList();  

    // Access to diffusion class. If we have no more species, we also can't diffuse
    getDiffusion()->mf_removeFromList(m_site);
  }
  
  m_site->m_updateNeighbours();

  // Remove count from list of sites that have n number of neighbours
  updateSiteCounter(numNeighbours, false);
  
}

void Desorption::mf_removeFromList() 
{ 
  m_lDesSites.remove(m_site); 
  //TODO: Is this necessary?
  m_site->removeProcess( this ); 
}

void Desorption::mf_addToList(Site *s) 
{ 
  m_lDesSites.push_back(s); 
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

  else if (m_site->getSpeciesName().size() == 0)
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

Diffusion* Desorption::getDiffusion()
{
  return m_pDiffusion;
}

void Desorption::setDiffusionPointer(Diffusion* d)
{
  m_pDiffusion = d;
}

bool Desorption::canDiffuse()
{
  return m_canDiffuse;
}

void Desorption::setDiffusion(bool canDiffuse)
{
  m_canDiffuse = canDiffuse;
}

void Desorption::updateSiteCounter(int neighbours, bool addOrRemove)
{
  //TODO: better bound checking  
  // Updates list of number of neighbours each possible site has
  if (addOrRemove)
  {
    m_numNeighbours.at(neighbours-1)++;
  }
  else
  {
    if (m_numNeighbours[neighbours-1] > 0)
    {
      m_numNeighbours[neighbours-1]--;
    }
  }
}

void Desorption::updateNeighbours(Site* s)
{
  // For all the neighbours of this site remove num sites and update m_numNeighbours
  vector<Site*> sites = s->getNeighs();
  //TODO: ensure getNeighs is properly updated within site
  for(vector<Site*> :: iterator itr = sites.begin(); itr != sites.end(); ++itr)
  {
    Site* site = *itr;
    updateSiteCounter(site->getNeighboursNum(), false);
    m_site = site;
    // Recalculate the number of sites
    s->m_updateNeighbours();
    updateSiteCounter(m_site->getNeighboursNum(), true);
  }
}

void Desorption::setSite(Site* s)
{
  m_site = s;
}

}
