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

#ifndef BCC_H
#define BCC_H

#include <iostream>
#include <stdlib.h>
#include <map>
#include <list>
#include <fstream>

#include "lattice.h"

using namespace std;
using namespace SurfaceTiles;
using namespace Utils;

class BCC : public Lattice
{
public:
  /// Constructor
  BCC(Apothesis *apothesis);

  /// Constructor
  BCC(Apothesis *apothesis, bool step, vector<int> stepInfo);

  /// Distructor.
  virtual ~BCC();

  void setSteps(bool hasSteps);

  void setStepInfo(int sizeX, int sizeY, int sizeZ);

  void mf_buildSteps();

  /// Sets the type of the lattice.
  void setType(string);

  /// Returns the x dimension of the lattice.
  inline int getX() { return m_iSizeX; }

  /// Returns the y dimension of the lattice.
  inline int getY() { return m_iSizeY; }

  /// Returns the size of the lattice.
  inline int getSize() { return m_iSizeX * m_iSizeY; }

  /// Returns a site with a specific id.
  Site *getSite(int id);

  /// Various checks if the lattice has been constucted correctly. Partially implemented.
  void check();

  /// Init the lattice.
  void init();

  /// Build the lattice with an intitial height.
  void build();

  /// Sets the minimun initial height for the lattice.
  void setInitialHeight(int height);

  /// Update neighbours
  void updateNeighbours(Site *s);

protected:
  /// Build the neighbours for the BCC lattice.
  void mf_neigh();

  /// Build the neighbours of each site depending on the type of the.
  void mf_buildNeighbours();

private:
  bool m_hasSteps;

  vector<int> m_stepInfo;
};

#endif // LATTICE_H
