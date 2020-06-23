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
#include "register.cpp"
#include "parameters.h"
#include <cmath>

namespace MicroProcesses{

/// Constructor
Diffusion::Diffusion
(
    vector<string> species,
    vector<double> energy,
    vector<double> frequency
)
:
m_sName("Diffusion"),
m_iNeighNum(0), 
m_apothesis(0),
m_diffusionSpecies(species),
m_diffusionEnergy(energy),
m_diffusionFrequency(frequency)
{
;  
}

Diffusion::~Diffusion(){;}

void Diffusion::init(){;}

//void Diffusion::setName( string name ){;}

string Diffusion::getName(){ return m_sName; }

void Diffusion::activeSites( Lattice* lattice){
  m_pLattice = lattice;
  vector< Site* > vSites = m_pLattice->getSites();

  for ( int i = 0; i < m_pLattice->getSize(); i++)
    if ( vSites[ i ]->getID()%2 != 0 ) {
      m_lAdsSites.push_back( vSites[ i ] );
      vSites[ i ]->addProcess( this );
    }
}

void Diffusion::selectSite()
{
  /* This comes from random i.e. picking from the available list for diffusion randomly */
  int y = rand()%getActiveList().size();
  int counter = 0;
  list<Site* >::iterator site = m_lAdsSites.begin();
  for (; site != m_lAdsSites.end(); site++ ){
    if ( counter == y )
      m_site = (*site);
    counter++;
    }
}

void Diffusion::perform()
{ 
  int height = m_site->getHeight();
  height = height + 2;
  m_site->setHeight( height);
  mf_removeFromList();
  mf_updateNeighNum();

  //this line is just to complete medium task for apothesis,
  //as it was written that the process should print it's name.
  //cout << "Hello, I am " << m_sName << endl;
}

void Diffusion::mf_removeFromList() { m_lAdsSites.remove( m_site); m_site->removeProcess( this ); }

void Diffusion::mf_addToList(Site *s) { m_lAdsSites.push_back( s); }

void Diffusion::mf_updateNeighNum()
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
double Diffusion::getProbability()
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

  double prob = -v0*exp((E-Em)/(dkBoltz*dTemp))*exp(-n*E/(dkBoltz*dTemp));
  return prob;
}

list<Site* > Diffusion::getActiveList()
{
  return m_lAdsSites;
}

void Diffusion::setProcessMap(map< Process*, list<Site* >* >* procMap )
  {
  m_pProcessMap = procMap;
  (*m_pProcessMap)[ this] = &m_lAdsSites;
  }

void Diffusion::test()
{
  cout << m_lAdsSites.size() << endl;
}

}
