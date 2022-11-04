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
#include "species_new.h"
#include <set>

using namespace std;
using namespace Utils;
using namespace SurfaceTiles;

class Lattice: public Pointers
  {
  public:
    /// The type of the lattice.
       enum Type{
               NONE,
               SimpleCubic,
               FCC
               };

    /// Constructor
    Lattice(Apothesis* apothesis);

    /// Distructor.
    virtual ~Lattice();

    /// Sets the type of the lattice.
    void setType( string );

    /// Sets the species on lattice sites.
    void setSpecies(Species*);

    /// Returns the x dimension of the lattice.
    inline int getX() { return m_iSizeX;} // = 0;

    /// Returns the y dimension of the lattice.
    inline int getY(){ return m_iSizeY;}// = 0;

    /// Returns the size of the lattice.
    inline int getSize(){ return m_iSizeX*m_iSizeY; } // = 0;

    /// Call update neighbours function;
//    virtual void updateNeighbours(Site* site) = 0;

//    /// Builds a  stepped surface
//    virtual void buildSteps(int iStepX, int iStepY);

    /// Returns the type of the lattice
    Lattice::Type getType();

    ///Returns the lattice
    Lattice* getLattice();

    /// Returns a site with a specific id.
    Site* getSite( int id);

    /// Returns a site with a specific id as in 2D space.
    Site* getSite( int i, int j);

    /// Returns all the sites of the lattice.
    vector<Site*> getSites();

    /// Init the lattice.
    void init();

    /// Set the type of the lattice.
    void setType( Type type );

    /// Set the X dimension of the lattice.
    void setX( int x );

    /// Set the Y dimension of the lattice.
    void setY( int y );

    /// Build the lattice with an intitial height.
    virtual void build() = 0;

    /// Sets the minimun initial height for the lattice.
    void setInitialHeight( int  height );

    //Set true if the lattice has steps
    inline void setSteps(bool hasSteps){m_hasSteps = hasSteps; }

    //Set true if the lattice has steps
    inline bool hasSteps(){ return m_hasSteps; }

    /// Set the "cut" of the surface
    void setOrientation(string s){ m_sOrient = s; }

    /// Store the surface step info
    void setStepInfo(int, int, int);

    //Prints the lattice heights
    void print();
    void printNeighNum();

    //Prints the neighbors
    void printNeighs(int);

    /// Returns the differnce between the first and last step
    inline int getStepDiff(){ return m_iStepDiff; }

    /// Write the lattice in XYZ format in a filename
    virtual void writeXYZ( string filename );

    virtual void writeLatticeHeights( double, int );

    /// Sets the steps of the surface in X
    inline void setStepX( int stepX ){ m_iStepX = stepX; }

    /// Sets the steps of the surface in Y
    inline void setStepY( int stepY ){ m_iStepY = stepY; }

    /// Returns the step in X
    inline int getStepX(){ return m_iStepX; }

    /// Returns the step in Y
    inline int getStepY(){ return m_iStepY; }

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

    /// True if the lattice has steps (comes from the input file if the Step keyword is found).
    bool m_hasSteps = false;

    //The info for the step surface i.e. [1 20 10]
    int m_iStepX;
    int m_iStepY;
    int m_iStepZ;

    map< string, set<int > >* m_pProcMap;

    Species * m_pSpecies;

    /// The cut of the lattice [ 100, 110, 111 ]
    string m_sOrient;

    /// The height differences between the first and last step
    int m_iStepDiff;
  };

#endif // LATTICE_H
