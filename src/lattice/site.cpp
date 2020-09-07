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

#ifndef SITE_CPP
#define SITE_CPP

#include "site.h"

namespace SurfaceTiles
{

Site::Site():
m_phantom(false)
{;}

Site::~Site(){;}

vector<Site* > Site::getNeighs() 
{ 
  return m_vNeigh; 
}

void Site::setNeigh( Site* s)
{
  m_vNeigh.push_back( s);
}

void Site::setID(int id) 
{ 
  m_iID = id;
}

int Site::getID() 
{ 
  return m_iID;
}

void Site::setHeight(int h)
{ 
  m_iHeight = h;
}

int Site::getHeight() 
{ 
  return m_iHeight; 
}

void Site::setNeighboursNum(int num)
{ 
  m_iNumNeighs = num;
}

int Site::getNeighboursNum()
{
  return m_vNeigh.size();
}
  
void Site::setNeighPosition(Site* s, NeighPoisition np)
{ 
  m_mapNeigh[np] = s; 
}

Site* Site::getNeighPosition(NeighPoisition np)
{ 
  return m_mapNeigh[np];
}

void Site::storeActivationSite(Site* s, ActivationSite as)
{
  m_mapAct[as] = s;
}

Site* Site::getActivationSite(ActivationSite as)
{ 
  return m_mapAct[as];
}

void Site::addSpecies(Species* s)
{
  m_species.push_back(s);
}

void Site::removeSpecies(Species* s)
{
  // TODO: Is this better to do by species pointer or by string?. I think pointer is easier
  // Search through species list until we've found the specific one to remove

  // Store previous state to ensure that we remove the desired species
  int numIter = 0;
  int m_species_prevSize = m_species.size();

  for (vector<Species*>::iterator itr = m_species.begin(); itr != m_species.end(); ++itr)
  {
    if (*itr == s)
    {
        m_species.erase(itr);
        break;
    }
    ++numIter;
  }
  
  // Output warning message if we didn't remove anything
  if (numIter == m_species_prevSize)
  {
    cout<<"Warning: did not find an instance of " << s->getName() << "in site " << getID()<<endl;
  }
}

vector<Species*> Site::getSpecies()
{
  return m_species;
}

vector<string> Site::getSpeciesName()
{
  vector<string> names;
  for (vector<Species*>::iterator itr = m_species.begin(); itr != m_species.end(); ++itr)
  {
    names.push_back((*itr)->getName());
  }
  return names;
}

void Site::addProcess(Process* process) 
{ 
  m_lProcs.push_back(process); 
}

void Site::removeProcess(Process* process) 
{
  m_lProcs.remove(process); 
}

list<Process*> Site::getProcesses()
{
  return m_lProcs;
}

void Site::setPhantom(bool phantom)
{
  m_phantom = phantom;
}

bool Site::isPhantom()
{
  return m_phantom;
}

void Site::removeDuplicates()
{
 //TODO
}
}


#endif
