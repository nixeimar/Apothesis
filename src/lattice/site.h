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
    void setHeight( int );

    /// Get the height of the particular site.
    int getHeight();

    /// Set the neigbours.
    void setNeigh( Site*);

    /// Initialize map
    void initSpeciesMap(int numSpecies);

    /// Get the neigbours at the same level.
    vector< Site* > getNeighs();

    /// Check if this site is active and can absorb.
    bool isActive();

    /// Set an ID for this site.
    void setID( int );

    /// Get the ID of this lattice site.
    int getID();

    /// Set the number of the neighbours according to the height of its neighbour sites
    inline void setNeighboursNum( int n ) { m_iNumNeighs = n; }

    /// Returns the number of the neighbours according to the height of its neighbour sites.
    int getNeighboursNum();

    /// Set the neihbour position for this site.
    void setNeighPosition( Site*, NeighPoisition  );

    /// Get the neihbour position for this site.
    Site* getNeighPosition( NeighPoisition );

    /// Stores the sites that are activated.
    void storeActivationSite( Site* s, ActivationSite as);

    /// Returns the sites that are activated.
    Site* getActivationSite( ActivationSite as);

    /// This will holds all the elements that can interact with the surfaces.
    /// Important when we talk about surface reactions.
    /// Element class has not implement yet for that we use forward decleration.
    void addSpecies( Species* s);

    void removeSpecies( Species* s);

    vector<Species*> getSpecies();

    vector<string> getSpeciesName();

    /// Add a processes in the list of processes that this site can participate in.
    void addProcess( Process* );

    /// Remove a processes from the list of processes that this site can participate in.
    void removeProcess( Process* );

    /// Get pointer to possible processes that can occur on this site
    list< Process* > getProcesses();

    // Set the lattice type that this site belongs to
    //void setLatticeType(LatticeType);

    /// Set the site as phantom (or not)
    void setPhantom(bool phantom);

    /// Return boolean of current state of site
    bool isPhantom();

    /// Return pointer to map
    map<int, int> getSpeciesMap();

    // Update neighbour list
    void m_updateNeighbours();

    // Call updateneighbours on another site
    void m_updateNeighbourList();

    // Reserve memory for the vector holding the species indexes for this site
    void initSpecies( int );

    // Get the react species
    inline valarray< int > getReactSpecies() { return m_vReacSpecies; }

    //Increase the height of the site by one
    void increaseHeight( int );

    //Decrease the height of the site by one
    void decreaseHeight( int );

    //Increase the number of neighs
    void increaseNeighsNum();

    /// Set the first negihbors of this site
    void set1stNeibors( int level, Site* );

    /// Returns the 1st neigbors
    inline map<int, vector<Site* > > get1stNeihbors() const { return m_m1stNeighs; }

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
};

}

#endif
