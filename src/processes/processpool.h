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

#ifndef PRORCESSPOOL_H
#define PRORCESSPOOL_H

#include <string>
#include <map>

#include "process.h"

using namespace std;

namespace MicroProcesses
{

class ProrcessPool
{
public:
    ProrcessPool();

    inline Process* getProcessByName(string name) { return m_mapProcs[ name ]; }
    inline Process* getProcessByID( int id ){ return m_mapProcsIDs[ id ]; }

    void addProcess( string name, Process* proc);
    void addProcess( int id, Process* proc);
    inline int getProcessNum() { return m_mapProcs.size(); }

private:
    map<string, Process*> m_mapProcs;
    map<int, Process*> m_mapProcsIDs;
};

}

#endif // PRORCESSPOOL_H
