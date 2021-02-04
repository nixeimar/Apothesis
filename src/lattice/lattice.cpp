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

#include "lattice.h"
#include "read.h"

Lattice::Lattice(Apothesis *apothesis) : Pointers(apothesis)
{
  //Document input =
}

void Lattice::setType(string sType)
{
  if (sType == "FCC")
    m_Type = FCC;
  else if (sType == "BCC")
    m_Type = BCC;
  else
    m_Type = NONE;
}

void Lattice::setX(int x) { m_iSizeX = x; }

void Lattice::setY(int y) { m_iSizeY = y; }

void Lattice::setInitialHeight(int height) { m_iHeight = height; }


Lattice::~Lattice()
{
}

vector<Site *> Lattice::getSites()
{
  return m_vSites;
}

Lattice::Type Lattice::getType()
{
  switch (m_Type)
  {
  case FCC:
    return FCC;
  case BCC:
    return BCC;
  default:
    return NONE;
  }
}

Site *Lattice::getSite(int id) { return m_vSites[id]; }

double Lattice::mf_roughness()
{
  double roughness = 0;

  for (auto& site:m_vSites)
  {
    for (auto& neighbour: site->getNeighs())
    {
      roughness += abs(site->getHeight() - neighbour->getHeight());
    }
  }

  return 1 + roughness / (2*m_vSites.size());
}

double Lattice::getRoughness()
{
  return mf_roughness();
}

void Lattice::check()
{
  int k = 0;

  cout << "Checking lattice..." << endl;

  int test = 2;
  cout << test << ": ";
  cout << "W:" << getSite(test)->getNeighPosition(Site::WEST)->getID() << " ";
  cout << "Wu:" << getSite(test)->getNeighPosition(Site::WEST_UP)->getID() << " ";
  cout << "WD:" << getSite(test)->getNeighPosition(Site::WEST_DOWN)->getID() << " ";
  cout << "E:" << getSite(test)->getNeighPosition(Site::EAST)->getID() << " ";
  cout << "EU:" << getSite(test)->getNeighPosition(Site::EAST_UP)->getID() << " ";
  cout << "ED:" << getSite(test)->getNeighPosition(Site::EAST_DOWN)->getID() << " ";
  cout << "N:" << getSite(test)->getNeighPosition(Site::NORTH)->getID() << " ";
  cout << "S:" << getSite(test)->getNeighPosition(Site::SOUTH)->getID() << endl;

  cout << "Activation: " << endl;

  cout << "N:" << getSite(test)->getActivationSite(Site::ACTV_NORTH)->getID() << " ";
  cout << "S:" << getSite(test)->getActivationSite(Site::ACTV_SOUTH)->getID() << " ";
  cout << "E:" << getSite(test)->getActivationSite(Site::ACTV_EAST)->getID() << " ";
  cout << "W:" << getSite(test)->getActivationSite(Site::ACTV_WEST)->getID() << endl;
}
