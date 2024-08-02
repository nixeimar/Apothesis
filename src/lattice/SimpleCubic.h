//============================================================================
//    Apothesis: A kinetic Monte Calro (KMC) code for deposition processes.
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

#ifndef SIMPLECUBIC_H
#define SIMPLECUBIC_H

#include <iostream>
#include <stdlib.h>
#include <map>
#include <list>
#include <fstream>

#include "io.h"
#include "lattice.h"

using namespace std;
using namespace SurfaceTiles;
using namespace Utils;

class SimpleCubic : public Lattice
{
public:
  /// Constructor
  SimpleCubic();

  /// Distructor.
  virtual ~SimpleCubic();

  void setSteps(bool hasSteps);

  /// Sets the type of the lattice.
  void setType(string);

  /// Various checks if the lattice has been constucted correctly. Partially implemented.
  void check();

  /// Init the lattice.
  void init();

  /// Build the lattice with an intitial height.
  void build() override;

  /// Create stepped surface
  void buildSteps() override;

  /// Sets the minimun initial height for the lattice.
  void setInitialHeight(int height);

  /// Calculate the number of neighbor based on the height
  int calculateNeighNum( int id );

  void writeLatticeHeights( double, int );

  inline int getNumFirstNeihgs() override { return 5; }

  inline bool isStepped(){return m_bHasSteps;}

  unordered_map<string, double> computeCoverages( vector<string> species) override;

protected:
  /// Build the neighbours for the BCC lattice for each site.
  void mf_neigh();

  /// Build the neighbours of each site depending on the type of the.
  void mf_buildNeighbours();

private:

  Lattice* m_lattice;   

  bool m_bHasSteps = false;

  vector<int> m_stepInfo;

  int m_iMinNeigs;

  int m_iSiteNeighsNum;

  REGISTER_LATTICE( SimpleCubic )
};

#endif // LATTICE_H
