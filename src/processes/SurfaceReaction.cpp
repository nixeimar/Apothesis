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
	vector<double> stoichiometry,
	double energy,
	double preExpFactor
)
:
m_sName("SurfaceReaction"),
m_iNeighNum(0), 
m_stoichiometry(stoichiometry),
m_energy(energy),
m_preExpFactor(preExpFactor)
{
  // Variable to check if the input file is configured properly
  bool readPositive = false;
  if (m_stoichiometry[0] >= 0)
  {
    m_apothesis->pErrorHandler->error_simple_msg("First stoichiometric input is positive or 0.");
  }

  double previousValue = m_stoichiometry[0];
  int counter = 0;
  for (vector<double> :: iterator itr = m_stoichiometry.begin(); itr != m_stoichiometry.end(); ++itr)
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
      m_lAdsSites.push_back( vSites[ i ] );
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
}

void SurfaceReaction::mf_removeFromList() { m_lAdsSites.remove( m_site); m_site->removeProcess( this ); }

void SurfaceReaction::mf_addToList(Site *s) { m_lAdsSites.push_back( s); }

void SurfaceReaction::mf_updateNeighNum()
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

//this process is not complete.
double SurfaceReaction::getProbability()
{
  /* These are parameters values (I/O) */
  double dNavogadro = m_apothesis->pParameters->dAvogadroNum;
  //double dPres = m_apothesis->pParameters->getPressure();
  double dTemp = m_apothesis->pParameters->getTemperature();
  double dkBoltz = m_apothesis->pParameters->dkBoltz;

  //these parameters are now given constant values for now.
  //later they will be taken from the input file.
  double v0 = 0.5;
  double E = 27500*4.2/dNavogadro;
  double Em = 0;
  double n = 5;

  double prob = 1;
  return prob;
}

list<Site* > SurfaceReaction::getActiveList()
{
  return m_lAdsSites;
}

void SurfaceReaction::setProcessMap(map< Process*, list<Site* >* >* procMap )
  {
  m_pProcessMap = procMap;
  (*m_pProcessMap)[ this] = &m_lAdsSites;
  }

void SurfaceReaction::test()
{
  cout << m_lAdsSites.size() << endl;
}

}
