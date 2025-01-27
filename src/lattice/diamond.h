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

#ifndef DIAMOND_H
#define DIAMOND_H

#include <iostream>
#include <stdlib.h>
#include <map>
#include <list>
#include <fstream>

#include "lattice.h"

using namespace std;
using namespace SurfaceTiles;
using namespace Utils;

class Diamond:public Lattice
{
public:
    /// The type of the lattice.

    /// Constructor
    Diamond( Apothesis* apothesis );

    /// Distructor.
    virtual ~Diamond(){}

    /// Build the lattice with an intitial height.
    void build();

    /// Sets the minimun initial height for the lattice.
    void setInitialHeight( int  );

    void readHeightsFromFile() override;

    void readSpeciesFromFile() override;

    /// Builds a  stepped surface
    virtual void buildSteps();

    /// Write the lattice in XYZ format in a filename
    void writeXYZ( string );

    unordered_map<string, double > computeCoverages( vector<string> species ) override;

protected:
    /// Build the first neighbours for the diamond lattice.
    void mf_neigh();
};

#endif // DIAMOND
