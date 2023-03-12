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

#ifndef PROCESS_H
#define PROCESS_H

#include <iostream>
#include <string>
#include <map>
#include <any>
#include "lattice.h"
#include "site.h"
#include "extLibs/random_generator.h"
#include "parameters.h"
#include "errorhandler.h"

#include "factory_process.h"

using namespace std;
using namespace SurfaceTiles;
using namespace Utils;

/** The pure virtual class from which every other process is generated.*/
namespace MicroProcesses
{

class Process
{

public:
    Process();
    virtual ~Process();

    ///Get probability
    virtual double getRateConstant() = 0;

    /// Perform this process in the site and compute/store the affected sites
    virtual void perform( Site* ) = 0;

    /// The rules for this type of process e.g. the neighbour of site Site.
    virtual bool rules( Site* ) = 0;

    /// Initialization for this process (e.g. temperature, pressure, mole fraction etc.)
    /// This must be for every process according to the process
    virtual void init( vector<string> params ){ m_vParams = params; }

    /// Returns the sites that are affected by this process including the site that this process is performed.
    inline set<Site*> getAffectedSites() { return m_seAffectedSites; }

    inline void setName( string procName ){ m_sProcName = procName; }
    inline string getName(){ return  m_sProcName; }

    inline void setID( int id ){ m_iID = id; }
    inline int getID(){ return m_iID; }

    inline void setLattice( Lattice* lattice ){ m_pLattice = lattice; }

    /// Counts how many times this process happens
    inline void eventHappened(){ m_iHappened++; }

    /// Returns how many times this process happened
    int getNumEventHappened(){ return m_iHappened; }

    /// Set the random generator
    inline void setRandomGen( RandomGen::RandomGenerator* randgen ) { m_pRandomGen = randgen; }

    inline void setSysParams( Utils::Parameters* p) { m_pUtilParams = p; }
    inline void setErrorHandler( ErrorHandler* error ) { m_error = error; }

    inline void setUncoAccepted( bool isUncoAccepted) { m_bUncoAccept = isUncoAccepted; }
    inline bool isUncoAccepted() { return m_bUncoAccept; }


    inline void setNumNeighs( int i){ m_iNumNeighs = i;}
    inline int getNumNeighs(){ return m_iNumNeighs;}

    inline void setNumVacantSites( int i){ m_iNumVacant = i;}
    inline int getNumVacantSites(){ return m_iNumVacant;}

protected:

    ///Pointer to the lattice of the process
    Lattice* m_pLattice;

    ///The parameters of the system and constant values
    Utils::Parameters* m_pUtilParams;

    /// Error handler for the processes
    ErrorHandler* m_error;

    /// Vector storing the variables for this processs.
    /// The first position in the vector is always a string declaring the type (e.g. simple, arrhenius etc.)
    /// followed by the parameters needed for this process to perform
    vector<string> m_vParams;

    ///A list holding all affected sites from this process
    set<Site*> m_seAffectedSites;

    ///The random generator
    RandomGen::RandomGenerator* m_pRandomGen;

    ///The type of the process
    string m_sType;

    ///The probability value
    double m_dProb;

    /// The name of this prcess
    string m_sProcName;

    ///Set true if it is always possible
    bool m_bUncoAccept;

    /// Checks if the specific species is part of the growing film
    bool isPartOfGrowth( string name);

    /// The number of sites occupied by the process (default 1)
    int m_iNumSites;

    /// The number of neighbours of this process (default 1)
    int m_iNumNeighs;

    /// The number of vacant sites of this process (default 1)
    int m_iNumVacant;

private:
    /// The id of the process
    int m_iID;

    /// Counts the times that this processes happened
    int m_iHappened;
};
}

#endif // Process_H

