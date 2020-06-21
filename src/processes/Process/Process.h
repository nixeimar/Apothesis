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

#ifndef TEMP_PROCESS_H

#define TEMP_PROCESS_H

#include <string>
#include <vector>
#include <list>
#include <iostream>
#include <map>

#include "lattice.h"
#include "factory_process.h"
#include "site.h"

using namespace std;
using namespace SurfaceTiles;

/** The pure virtual class from which every other process is generated.*/

namespace MicroProcesses{

template<class ProcessType>
class Process;

template<class ProcessType>
class Process:
    public ProcessType
{
  public:
    /// Constructor of the interface.
    Process(){}

    /// Destructor.
    virtual ~Process(){}

    /// Set the name of the process.
    virtual void setName( string s) =0;

    /// Get the name of the process.
    virtual string getName() =0;

    /// Constructs the sites that a process can be performed
    virtual void activeSites( Lattice* ) =0;

    /// The site that this process will be performed.
    /// The site is selected from the available sites that have been constructed in activeSites
    virtual void selectSite() = 0;

    /// The process map which holds all the processes and the sites that each can be performed.
    // Not to handy. Re-think...
    virtual void setProcessMap( map< Process*, list<Site* >* >* ) =0;

    /// Perform the process.
    virtual void perform() = 0;

    /// Calculate and get the Probability of this process.
    virtual double getProbability() = 0;

    /// Get the list of active sites where the process can be performed.
    /// This is updated after a process is performed.
    virtual list<Site* > getActiveList() =0;

    /// Set the instance of kmc that this process will be performed.
    virtual void setInstance( Apothesis* apothesis ) = 0;
  };
}
};
#endif // PROCESS_H
