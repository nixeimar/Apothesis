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
    virtual void perform( int siteID ) = 0;

    /// The rules for this type of process e.g. the neighbour of site Site.
    virtual bool rules( Site* ) = 0;

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

protected:
    /** Pointer to the lattice of the process */
    Lattice* m_pLattice;

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

