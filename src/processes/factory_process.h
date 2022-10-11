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

#ifndef FACTORY_PROCESS_H
#define FACTORY_PROCESS_H

#include <string>
#include <map>

#include "register.h"

/** A factory process (pattern) for creating the processes
 * for the surface processes. */
namespace MicroProcesses {
class Process;
}

class AbstractProcess;

using namespace std;

class FactoryProcess
{
public:
    /// Constructor
    FactoryProcess();

    /// Destructor
    virtual ~FactoryProcess();

    /// Register the new process (Store it in the map table since a map table cannot be initialized)
    static void registerThis( const string& , AbstractProcess* );

    /// Create the process
    static MicroProcesses::Process* createProcess( const string& );

private:
    /// The factory map
    static map< string, AbstractProcess* >& getTable();

};

#define REGISTER_PROCESS( __NAME__ ) \
    private: static const Register< __NAME__> creator;

#define REGISTER_PROCESS_IMPL( __NAME__) \
    const Register< __NAME__ > __NAME__::creator(#__NAME__);

#endif // FACTORY_PROCESS_H
