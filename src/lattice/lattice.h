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

#ifndef LATTICE_H
#define LATTICE_H

#include <iostream>
#include <stdlib.h>
#include <map>
#include <list>
#include <fstream>

#include "pointers.h"
#include "site.h"
#include "errorhandler.h"

using namespace std;
using namespace SurfaceTiles;
using namespace Utils;

class Lattice: public Pointers
  {
  public:
    /// The type of the lattice.
       enum Type{ NONE,
               BCC,
               FCC
               };

    /// Constructor
    Lattice( Apothesis* apothesis );

    /// Distructor.
    virtual ~Lattice();

    /// Sets the type of the lattice.
    void setType( string );

    /// Returns the x dimension of the lattice.
    inline int getX(){return m_iSizeX;}

    /// Returns the y dimension of the lattice.
    inline int getY(){return m_iSizeY;}

    /// Returns the size of the lattice.
    inline int getSize(){return  m_iSizeX*m_iSizeY;}

    /// Returns the type of the lattice
    Lattice::Type getType();

    ///Returns the lattice
    Site* getLattice();

    /// Returns a site with a specific id.
    Site* getSite( int id);

    /// Returns all the sites of the lattice.
    vector<Site*> getSites();

    /// Various checks if the lattice has been constucted correctly. Partially implemented.
    void check();

    /// Init the lattice.
    void init();

    /// Set the type of the lattice.
    void setType( Type type );

    /// Set the X dimension of the lattice.
    void setX( int x );

    /// Set the Y dimension of the lattice.
    void setY( int y );

    /// Build the lattice with an intitial height.
    void build();

    /// Sets the minimun initial height for the lattice.
    void setInitialHeight( int  height );

    /// Returns random site
    Site* randomSite();

  protected:
    /// The size of the lattice in the x-dimension.
    int m_iSizeX;

    /// The size of the lattice in the y-dimension.
    int m_iSizeY;

    /// The minimum initialize size of the lattice.
    int m_iHeight;

    /// The type of the lattice: BCC, FCC etc.
    Type m_Type;

    /// The sites that consist the lattice.
    vector<Site* > m_vSites;

    /// The neighbours for the FCC lattice.
    void mf_fccNeigh();

    /// Build the neighbours of each site depending on the type of the.
    void mf_buildNeighbours();
  };

#endif // LATTICE_H
