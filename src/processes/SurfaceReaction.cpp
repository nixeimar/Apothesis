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

#include "SurfaceReaction.h"
#include "parameters.h"
#include "register.cpp"

namespace MicroProcesses{


/// Constructor
SurfaceReaction::SurfaceReaction
(
  Apothesis* apothesis,
  vector<Species*> species,
	vector<double> const& stoichiometry,
	double energy,
	double preExpFactor,
  bool immobilized
)
:
m_sName("SurfaceReaction"),
m_iNeighNum(0), 
m_stoichiometry(stoichiometry),
m_energy(energy),
m_preExpFactor(preExpFactor),
m_immobilized(immobilized)
{
  // Variable to check if the input file is configured properly
  bool readPositive = false;

  if (m_stoichiometry[0] >= 0)
  {
    m_apothesis->pErrorHandler->error_simple_msg("First stoichiometric input is positive or 0.");
  }
  double previousValue = m_stoichiometry[0];
  int counter = 0;
  for (vector<double> :: const_iterator itr = m_stoichiometry.begin(); itr != m_stoichiometry.end(); ++itr)
  {
    // Check to make sure first digit is negative, and no negative numbers appear after a positive
    double stoichiometry = *itr;
    if (stoichiometry > 0)
    {
      readPositive = true;
    }

    if (stoichiometry < 0 && readPositive)
    {
      m_apothesis->pErrorHandler->error_simple_msg("Stoichiometric configuration in input file is incorrect.");
    }
    if (stoichiometry < 0)
    {
      m_stoichReactants.push_back(-1 * stoichiometry);
      m_reactants.push_back(species[counter]);
    }
    else if (stoichiometry > 0)
    {
      m_stoichProducts.push_back(stoichiometry);
      m_products.push_back(species[counter]);
    }
    else
    {
      m_apothesis->pErrorHandler->error_simple_msg("Warning, value of stoichiometric coefficient is 0");
    }
    counter++; 
  }

}

SurfaceReaction::~SurfaceReaction(){;}

void SurfaceReaction::init(){;}


string SurfaceReaction::getName(){ return m_sName; }

void SurfaceReaction::activeSites( Lattice* lattice){
  m_pLattice = lattice;
  vector< Site* > vSites = m_pLattice->getSites();

  for ( int i = 0; i < m_pLattice->getSize(); i++)
    if ( vSites[ i ]->getID()%2 != 0 ) {
      vSites[ i ]->addProcess( this );
    }
}

void SurfaceReaction::selectSite()
{
  /* This comes from random i.e. picking from the available list for SurfaceReaction randomly */
  int y = rand()%getActiveList().size();
  int counter = 0;
  list<Site* >::iterator site = m_lAdsSites.begin();
  for (; site != m_lAdsSites.end(); site++ ){
    if ( counter == y )
      m_site = (*site);
    counter++;
    }
}

void SurfaceReaction::perform()
{ 
  // Perform surface reaction here
  vector<Species*> :: iterator rItr = m_reactants.begin();
  vector<Species*> :: iterator pItr = m_products.begin();
  for (; rItr != m_reactants.end(); ++rItr)
  {
    for (int i = 0; i < (*rItr)->getStoicCoeff()*-1; ++i)
    {
      // TODO Call desorption, or simply remove?
      m_site->removeSpecies(*rItr);
    }
  }
  for (; pItr != m_products.end(); ++pItr)
  {
    for (int i = 0; i < (*pItr)->getStoicCoeff(); ++i)
    {
      // Add species to site. Need to do anything else?
      m_site->addSpecies(*pItr);
    }
    m_site->setHeight(m_site->getHeight()+2); //TODO generalize
  }

  if (!canReact(m_site))
  {
    m_activeSites--;
  }
  //TODO 
  //if (m_immobilized)
  //{
  //  vector<Adsorption*> pAds = m_apothesis->getAdsorptionPointers();
  //  for (vector<Adsorption*> :: iterator itr = pAds.begin(); itr != pAds.end(); ++itr)
  //  {
  //    // Allow adsorption once more
  //    Adsorption* a = *itr;
  //    a->mf_addToList(m_site);
  //  }
  //}
}

void SurfaceReaction::mf_removeFromList() { m_lAdsSites.remove( m_site); m_site->removeProcess( this ); }

void SurfaceReaction::mf_addToList(Site *s) { m_lAdsSites.push_back( s); }


//this process is not complete.
double SurfaceReaction::getProbability()
{
  if (m_lAdsSites.size() < 1)
  {
    return 0;
  }

  Parameters* parameters = m_apothesis->pParameters;
  double rate = m_preExpFactor * exp(-m_energy/parameters->getTemperature()/parameters->dR);
  delete parameters;
  return rate * m_activeSites;
}

list<Site* > SurfaceReaction::getActiveList()
{
  return m_lAdsSites;
}

const vector<double> SurfaceReaction::getStoichiometry()
{
  return m_stoichiometry;
}

void SurfaceReaction::test()
{
  cout << "size in test is " << m_stoichiometry.size()<<endl;
}

bool SurfaceReaction::canReact(Site* site)
{

  vector<Species*> species = site->getSpecies();
  vector<Species*> :: iterator sItr = species.begin();

  map<int, int> map_species = site->getSpeciesMap();
  
  bool canReact = false;

  // Needs double checking! How to avoid re-computation at every itr?
  for(int i = 0; i < m_reactants.size(); ++i)
  {
    if (map_species.at(m_reactants[i]->getId()) < m_stoichReactants[i])
    {
      return false;
    }
  }
  m_activeSites++;
  mf_addToList(site);
  cout<<"can react!"<<endl;
  return true;
}

void SurfaceReaction::setProcessMap(map< Process*, list<Site* >* >* procMap )
{
  
}

}
