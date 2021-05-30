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

#ifndef SITE_H
#define SITE_H

#include <string>
#include <vector>
#include <list>
#include <map>
#include <valarray>

#include "process.h"
#include "species.h"

using namespace std;
using namespace MicroProcesses;

/**  The site is where a process will be performed. The lattice is
 * a series of sites put together in space with certain symmetry. */

namespace SurfaceTiles
{

class Site
{
public:

    /// Contructor.
    Site();

    /// Destructor.
    virtual ~Site();

    //This should be on the Lattice
    //and every lattice build their own Neghbours
    /*    enum LatticeType {
                        BCC,
                        FCC
                     };*/

    /// The position of the neighbour.
    // This should be SOUTH_EAST, SOUTH_WEST etc ..
    enum NeighPoisition{
        WEST,
        WEST_UP,
        WEST_DOWN,
        EAST,
        EAST_UP,
        EAST_DOWN,
        NORTH,
        SOUTH
    };

    /// The site that is activated if the neighbours are all occupied.
    enum ActivationSite{
        ACTV_WEST,
        ACTV_WEST_UP,
        ACTV_WEST_DOWN,
        ACTV_EAST,
        ACTV_EAST_UP,
        ACTV_EAST_DOWN,
        ACTV_NORTH,
        ACTV_SOUTH
    };

    /// Set the height of the particular site.
    inline void setHeight( int h ) { m_iHeight = h; }

    /// Get the height of the particular site.
    inline int getHeight() { return m_iHeight; }

    /// Set the neigbours.
    inline void setNeigh(Site *s){ m_vNeigh.push_back(s); }

    /// Get the neigbours at the same level.
    inline vector<Site *> getNeighs() { return m_vNeigh; }

    /// Set an ID for this site.
    inline void setID(int id) { m_iID = id; }

    /// Get the ID of this lattice site.
    inline int getID() { return m_iID; }

    /// Set the number of the neighbours according to the height of its neighbour sites
    inline void setNeighsNum( int n ) { m_iNumNeighs = n; }

    /// Returns the number of the neighbours according to the height of its neighbour sites.
    inline int getNeighsNum(){ return m_iNumNeighs; }

    /// Set the neihbour position for this site.
    inline void setNeighPosition(Site *s, NeighPoisition np) { m_mapNeigh[ np ] = s; }

    /// Get the neihbour position for this site.
    Site* getNeighPosition(NeighPoisition np){ return m_mapNeigh[np]; }

    /// This will holds all the elements that can interact with the surfaces.
    /// Important when we talk about surface reactions.
    void addSpecies( Species* s);

    /// Remove a species from the site
    void removeSpecies( Species* s);

    /// Get the species currently in the site
    vector<Species*> getSpecies();

    vector<string> getSpeciesName();

    // Get the react species
    inline valarray< int > getReactSpecies() { return m_vReacSpecies; }

    /// Increase the height of the site by one
    inline void increaseHeight( int i ){ m_iHeight += i; }

    /// Decrease the height of the site by one
    inline void decreaseHeight( int i ){ m_iHeight -= i; }

    /// Set the first negihbors of this site
    void set1stNeibors( int level, Site* s) { m_m1stNeighs.at( level ).push_back( s ); }

    /// Returns the 1st neigbors
    inline map<int, vector<Site* > > get1stNeihbors() const { return m_m1stNeighs; }

    /// Returns true if is in lower step (used in the step case only)
    void setLowerStep( bool b){ m_isLowerStep = b; }

    /// Returns true if is in higher step (used in the step case only)
    void setHigherStep( bool b) { m_isHigherStep = b; }

    /// Returns true if is in lower step (used in the step case only)
    bool isLowerStep(){ return m_isLowerStep; }

    /// Returns true if is in higher step (used in the step case only)
    bool isHigherStep() { return m_isHigherStep; }

protected:
    //The lattice type that this site belongs to
    //LatticeType m_LatticeType;

    /// The ID of the site.
    int m_iID;

    /// The indexes coming from the vector.
    int m_iIIndex;
    int m_iJIndex;

    /// The height in the particular position.
    int m_iHeight;

    /// The neighbours at the same level - We do not need this because we have the map (see below).
    vector< Site*> m_vNeigh;

    /// Holds the number of the neighbours of the particular site according to each height/
    int m_iNumNeighs;

    //TODO: get rid of map?
    /// A map that holds the neighbour sites according to their orientation.
    map< NeighPoisition, Site*> m_mapNeigh;

    /// A map that hold the sites that this site can activate if it has all the its neighbours occupied.
    map< ActivationSite, Site* > m_mapAct;

    /// A list of active sites for each process
    map< Process*, vector<Site*>> activeSites;

    /// A map counting the species that are present
    map<int, int> m_mapSpecies;

    /// To check if a site is occupied or not.
    bool m_bIsOccupied;

    /// The list of processes that this site can participate in.
    list< Process* > m_lProcs;

    vector<Species* > m_species;

    bool m_phantom;

    //Holds the number of the reactants species in this site.
    // m_vReacSpecies has size the number of ALL reactants species
    valarray< int > m_vReacSpecies;

    /// The number of first neighs
    int m_iFirstNeighs;

    /// The number of second neighs
    int m_iSecondNeighs;

private:
    /// Orientation of the lattice
    string m_sOrint;

    /// The species that are in this site
    vector<species_new* > m_vSpecies;

    /// The 1st neighbors in the different levels
    /// below level
    /// same level
    /// upper level
    map< int, vector <Site* > > m_m1stNeighs;

    /// The 2nd nearest neighbors in the different levels
    /// below level
    /// same level
    /// upper level
    map< int, vector <Site* > > m_m2ndNeighs;

    /// if the site belong to the lower step storing vaious info
    bool m_isLowerStep;

    /// if the site belong to the higher step storing vaious info
    bool m_isHigherStep;
};

}

#endif
