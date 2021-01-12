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

  Site::Site(Lattice *lattice) : m_phantom(false),
                                 m_lattice(lattice)
  {
  }

  Site::~Site() { ; }

  vector<Site *> Site::getNeighs()
  {
    return m_vNeigh;
  }

  void Site::setNeigh(Site *s)
  {
    m_vNeigh.push_back(s);
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

  void Site::setLatticeType(LatticeType type)
  {
    m_LatticeType = type;
  }

  int Site::getNeighboursNum()
  {
    return m_vNeigh.size();
  }

  void Site::setNeighPosition(Site *s, NeighPoisition np)
  {
    m_mapNeigh[np] = s;
  }

  Site *Site::getNeighPosition(NeighPoisition np)
  {
    return m_mapNeigh[np];
  }

  void Site::storeActivationSite(Site *s, ActivationSite as)
  {
    m_mapAct[as] = s;
  }

  Site *Site::getActivationSite(ActivationSite as)
  {
    return m_mapAct[as];
  }

  void Site::addSpecies(Species *s)
  {
    m_species.push_back(s);
    m_mapSpecies[s->getId()]++;
  }

  void Site::removeSpecies(Species *s)
  {
    // TODO: Is this better to do by species pointer or by string?. I think pointer is easier
    // Search through species list until we've found the specific one to remove

    // Store previous state to ensure that we remove the desired species
    int numIter = 0;
    int m_species_prevSize = m_species.size();

    int numOfSpecies = m_mapSpecies[s->getId()];

    if (numOfSpecies > 0)
    {
      for (vector<Species *>::iterator itr = m_species.begin(); itr != m_species.end(); ++itr)
      {
        if (*itr == s)
        {
          m_species.erase(itr);
          break;
        }
        ++numIter;
      }
      // Decrement number of said species
      m_mapSpecies[s->getId()]--;
    }

    // Output warning message if we didn't remove anything
    if (numOfSpecies < 1)
    {
      cout << "Warning: did not find an instance of " << s->getName() << "in site " << getID() << endl;
    }
  }

  vector<Species *> Site::getSpecies()
  {
    return m_species;
  }

  vector<string> Site::getSpeciesName()
  {
    vector<string> names;
    for (vector<Species *>::iterator itr = m_species.begin(); itr != m_species.end(); ++itr)
    {
      names.push_back((*itr)->getName());
    }
    return names;
  }

  void Site::addProcess(Process *process)
  {
    m_lProcs.push_back(process);
  }

  void Site::removeProcess(Process *process)
  {
    m_lProcs.remove(process);
  }

  list<Process *> Site::getProcesses()
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

  map<int, int> Site::getSpeciesMap()
  {
    return m_mapSpecies;
  }

  void Site::m_updateNeighbours()
  {
    int siteHeight = getHeight();
    int totalNeigh = 0;
    m_vNeigh.clear();

    
    // Check NESW sites, see if the heights are the same. If same, add to list of neighbours.
    bool isActiveEAST = false;
    isActiveEAST = (siteHeight == getNeighPosition(Site::EAST)->getHeight() &&
                    siteHeight == getNeighPosition(Site::EAST_DOWN)->getHeight() &&
                    siteHeight == getNeighPosition(Site::EAST_UP)->getHeight());
    if (isActiveEAST)
    {
      m_vNeigh.push_back(getNeighPosition(Site::EAST));
      totalNeigh++;
    }

    bool isActiveWEST = false;
    isActiveWEST = (siteHeight == getNeighPosition(Site::WEST)->getHeight() &&
                    siteHeight == getNeighPosition(Site::WEST_DOWN)->getHeight() &&
                    siteHeight == getNeighPosition(Site::WEST_UP)->getHeight());
    if (isActiveWEST)
    {
      m_vNeigh.push_back(getNeighPosition(Site::WEST));
      totalNeigh++;
    }

    bool isActiveNORTH = false;
    isActiveNORTH = (siteHeight == getNeighPosition(Site::WEST_UP)->getHeight() &&
                     siteHeight == getNeighPosition(Site::EAST_UP)->getHeight() &&
                     siteHeight == getNeighPosition(Site::NORTH)->getHeight());
    if (isActiveNORTH)
    {
      m_vNeigh.push_back(getNeighPosition(Site::NORTH));
      totalNeigh++;
    }

    bool isActiveSOUTH = false;
    isActiveSOUTH = (siteHeight == getNeighPosition(Site::WEST_DOWN)->getHeight() &&
                     siteHeight == getNeighPosition(Site::EAST_DOWN)->getHeight() &&
                     siteHeight == getNeighPosition(Site::SOUTH)->getHeight());
    if (isActiveSOUTH)
    {
      m_vNeigh.push_back(getNeighPosition(Site::SOUTH));
      totalNeigh++;
    }
  }

  void Site::m_updateNeighbourList()
  {
    // If called on a different site, simply call the respective site's m_updateNeighbours() function
    vector<Site *> neighbours{getNeighPosition(Site::EAST), getNeighPosition(Site::WEST), getNeighPosition(Site::NORTH), getNeighPosition(Site::SOUTH)};
    for (int neigh = 0; neigh < neighbours.size(); ++neigh)
    {
      neighbours[neigh]->m_updateNeighbours();
    }
  }

  void Site::initSpeciesMap(int numSpecies)
  {
    for (int i = 0; i < numSpecies; ++i)
    {
      m_mapSpecies[i] = 0;
    }
  }

} // namespace SurfaceTiles

#endif
