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

#include "processpool.h"

namespace MicroProcesses
{

ProrcessPool::ProrcessPool(){}

void ProrcessPool::addProcess( string name, Process* proc)
{
    pair<string, Process* > p;
    p.first = name;
    p.second = proc;

    m_mapProcs.insert( p );
}

void ProrcessPool::addProcess( int id, Process* proc)
{
    pair<int, Process* > p;
    p.first = id;
    p.second = proc;

    m_mapProcsIDs.insert( p );
}




}
