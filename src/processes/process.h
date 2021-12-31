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

#ifndef PROCESS_H
#define PROCESS_H

#include <iostream>
#include <string>
#include <map>
#include <any>
#include "lattice.h"
#include "site.h"
#include "species_new.h"
#include "aux/random_generator.h"

#include "factory_process.h"

using namespace std;
using namespace SurfaceTiles;

/** The pure virtual class from which every other process is generated.*/
namespace MicroProcesses
{

class Process
{

public:
    Process();
    virtual ~Process();

    ///Get probability
    virtual double getProbability() = 0;

    /// Perform this process in the site
    virtual void perform( Site* ) = 0;

    /// The rules for this type of process e.g. the neighbour of site Site.
    virtual bool rules( Site* ) = 0;

    /// Initialization for this process (e.g. temperature, pressure, mole fraction etc.)
    /// This must be for every process
    void init( map<string, any> params ){ m_mParams = params; }

    /// Returns the sites that are affected by this process including the site that this process is performed.
    inline set<Site*> getAffectedSites() { return m_seAffectedSites; }

    inline void setName( string procName ){ m_sProcName = procName; }
    inline string getName(){ return  m_sProcName; }

    inline void setID( int id ){ m_iID = id; }
    inline int getID(){ return m_iID; }

    inline void setTargetSite( int id );

    inline void setLattice( Lattice* lattice ){ m_pLattice = lattice; }

    /// Counts how many times this process happens
    inline void eventHappened(){ m_iHappened++; }

    /// Returns how many times this process happens
    int getNumEventHappened(){ return m_iHappened; }

    /// Insert a parameter for this process
    inline void setParameter( string str, any val ) { m_mParams[ str ] = val; }

    /// Return a parameter of this process
    any getParameter( string str ) { return m_mParams[ str ]; }

    /// Set the random generator
    inline void setRandomGen( RandomGen::RandomGenerator* randgen ) { m_pRandomGen = randgen; }

    inline void setUncoAccepted( bool isUncoAccepted) { m_bUncoAccept = isUncoAccepted; }
    inline bool isUncoAccepted() { return m_bUncoAccept; }

    inline void setConstant( bool isConst =false ){ m_bConstant = isConst; }

    //If this process is supposed to be constant (i.e. independent of coverage)
    inline bool isConstant(){ return m_bConstant;}

    inline void setNeighs( int num ){ m_iNeighs = num; }
    inline int getNeighs(){ return m_iNeighs; }

protected:
    /** Pointer to the lattice of the process */
    Lattice* m_pLattice;

    /// Map for storing the variables for this processs
    map<string, any> m_mParams;

    ///A list holding all affected sites from this process
    set<Site*> m_seAffectedSites;

    RandomGen::RandomGenerator* m_pRandomGen;

    bool m_bUncoAccept;

    bool m_bConstant;

    int m_iNeighs;

private:
    /// The name of this prcess
    string m_sProcName;

    /// The id of the process
    int m_iID;

    ///The type of the process
    string m_sType;

    /// Counts the times that this processes happened
    int m_iHappened;
};
}

#endif // Process_H

