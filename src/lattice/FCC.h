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

#ifndef FCC_H
#define FCC_H

#include <iostream>
#include <stdlib.h>
#include <map>
#include <list>
#include <fstream>

#include "lattice.h"
// #include "register_lattice.h"


using namespace std;
using namespace SurfaceTiles;
using namespace Utils;

class FCC: public Lattice
{ 
public:
    /// The type of the lattice.

    /// Constructor
    FCC();

    /// Distructor.
    virtual ~FCC();     

    /// Sets the type of the lattice.
    // void setType( string );

    /// Various checks if the lattice has been constucted correctly. Partially implemented.
    void check();

    /// Init the lattice.
    void init();

    /// Build the lattice with an intitial height.
    void build();

    /// Sets the minimun initial height for the lattice.
    void setInitialHeight( int  );

    /// Count the neighbors in the different levels
    int calculateNeighNum( int,  int);

    /// Write the lattice in XYZ format in a filename
    void writeXYZ( string );

    void buildSteps(int, int);

    unordered_map<string, double > computeCoverages( vector<string> species ) override;

private:     
    Lattice* m_lattice;   

protected:
    /// Build the first neighbours for the FCC(100) lattice.
    void mf_neigh_100();

    /// Build the first neighbours for the FCC(110) lattice.
    void mf_neigh_110();

    /// Build the first neighbours for the FCC(111) lattice.
    void mf_neigh_111();

    REGISTER_LATTICE(FCC)
};

#endif // LATTICE_H
